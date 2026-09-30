%%%%%%%
% Code to simulate filament swarms

% Plotting tool obtained from Matlab file exchange, used under the following license:
% 
% Copyright (c) 2016, Janus H. Wesenberg
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


clear;


%define initial condition
inputIC = readlines('input_initial_condition.txt');
inputIC = inputIC(1);
load(inputIC, 'W');
W0 = W(:,:,end);
worms = size(W0,1)/3;
nw = size(W0,2)-2; %size of ribbons


%name of output file
fileLabel = strcat('results/output_sim');

%Kirchhoff equation parameters
inpElastic = readmatrix('input_simulation_parameters.txt', 'range', 'A1:A12');

%target curvature
prefK01 = inpElastic(1);

%number of timesteps for the preferred curvature to reach its terminal value
Nbend = inpElastic(2); %either 0 or 25000

%elastic parameters
h = inpElastic(3); %radius
Eb = inpElastic(4);  %bending modulus
rho = inpElastic(5); %density

%damping and friction coefficient
eta = inpElastic(6);
zeta = inpElastic(7);
skindrag = inpElastic(8);

%brownian parameters
driftForce = inpElastic(9); %Magnitude of random force

%contact parameters
Kbulk = inpElastic(10); %sets ribbon-ribbon contact forces
bdryForceStrength = inpElastic(11); %sets force from boundary
fStrength = inpElastic(12); %inverse friction parameter, appears in boundary contact force




%%derived parameters
%elastic
Ar=pi*h^2; %cross-sectional area
I=0.25*pi*h^4; %moment of inertia
alpha = Eb*I; %bending energy scale
gamma = Eb*Ar; %stretching energy scale

%damping
eta_eff = eta * Ar;  %second derivative damping
zeta_eff = zeta*Ar; %collision stick factor


%Timing parameters which can be changed
N=50000; %number of timesteps, 50000
framenumber=200; %number of frames to capture, 200
dt= 10^(-4); %size of timesteps


%preferred curvatures
prefK1 = [linspace(0,prefK01,Nbend), prefK01+zeros(1,N-Nbend)];
prefK2 = [zeros(1,N)];

%container dimensions (possibly time-dependent)
diskR = 10 + [zeros(1,N)]; %radius of container
tubeL = 20 + [zeros(1,N)]; %height of the container




%%%%%%%%%%%%%%%%%
%Initialize variable for simulation
%%%%%%%%%%%%%%%%%

%initialize geometric parameters
[L0w, Ls0w, vLengths_fric] = initialLengths(W0);
ds = vLengths_fric';
massList = rho*Ar*vLengths_fric(:)';

%initialize contacts
V0comp = zeros(3, worms*(nw+2));
for ii = 1:worms
    V0comp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = W0(3*(ii-1)+1:3*ii,:);
end
Rw = dists(V0comp);
for ii = 1:worms-1 %get rid of ghost links between worm i head and worm i+1 tail
    Rw(:,ii*(nw+2)) = 100*h; Rw(ii*(nw+2),:) = 100*h;
end
%effective radius
iMw0 = interactMat(Rw,1.25*h);

%material frame
d1w = initialMaterialFrame(W0);

%Initialize filament position array, orientation, force and drift direction
W = zeros(3*worms,nw+2,framenumber+1); W(:,:,1) = W0;
velw = zeros(3,size(W,2),worms); velW = zeros(3, size(W,2)*worms);
accelw = zeros(3,size(W,2),worms); accelW = zeros(3, worms*(nw+2));
d1wList = zeros(3*worms, nw+1, framenumber+1); d1wList(:,:,1) = d1w;

V0 = W0; 
Ftot = zeros(3, nw+2, worms);

% main loop
tic; LoopTime0 = 0;
for i=1:N 

    %initialize variables, timestepping takes V00 -> V0, etc.
    V00 = V0;
    for ii = 1:worms
        w00 = V00(3*ii-2:3*ii,:);

        V0(3*ii-2:3*ii,:) = V00(3*ii-2:3*ii,:) + dt*velw(:,:,ii) + 0.5*dt^2*accelw(:,:,ii);

        d1w(3*(ii-1)+1:3*ii,:) = fupdate(V00(3*(ii-1)+1:3*ii,:), V0(3*(ii-1)+1:3*ii,:), d1w(3*(ii-1)+1:3*ii,:)); %update reference frame along the rope
    end


    %contact handling terms
    V0comp = zeros(3, worms*(nw+2));
    for ii = 1:worms
        V0comp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = V0(3*(ii-1)+1:3*ii,:);
    end

    Rw = dists(V0comp);
    for ii = 1:worms-1 %get rid of links between ribbon i head and ribbon i+1 tail
        Rw(:,ii*(nw+2)) = 100*h; Rw(ii*(nw+2),:) = 100*h; %no self-effects for rigid filaments
        
    end


    %Nonlinear adhesion
    iMw1 = interactMat(Rw,1.25*h); %change in effective radius here
    iMw = (iMw0+iMw1)>0;
    pMatw = pressureThreshold(Rw,iMw,1.25*h,Kbulk,0); %change in effective radius here
    iMw0 = iMw1; %update interaction matrix for next run

    %calculate forces and torques
    Ftot = zeros(3, nw+2, worms);
    iFW = iForce(V0comp,pMatw,Rw);
    for ii=1:worms

        w  =  V0(3*(ii-1)+1:3*ii,:);
        d1 = d1w(3*(ii-1)+1:3*ii,:);
        
        %elastic forces
        Fbend = bForce_aniso1(w, d1, prefK1(i), alpha) + bForce_aniso2(w, d1, prefK2(i), alpha);
        Felast = Fbend + sForce(w,L0w(ii,:),gamma); 

        %contact forces
        iF = iFW(:, (ii-1)*(nw+2)+1: ii*(nw+2));
        Fcont = contactForce(w, diskR(i), tubeL(i), h, bdryForceStrength, fStrength);
        
        Ftot(:,:,ii) = iF + Fcont+ Felast + sqrt(2*driftForce/dt) * sqrt(ds(ii,:)).*normrnd(0,1,3,nw+2); %total force on ribbon

    end


    %Update velocity and acceleration
    Fcomp = zeros(3, worms*(nw+2));
    for ii = 1:worms
        Fcomp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = Ftot(:,:,ii);
    end
    [velW, accelW] = velUpdate(worms, dt, Fcomp, velW, accelW, V0comp, iMw, vLengths_fric, massList, L0w, eta_eff, zeta_eff, skindrag, rho );
    velw = zeros(3,nw+2, worms);
    accelw = zeros(3,nw+2, worms);
    for ii=1:worms
        velw(:,:,ii) = velW(:,(nw+2)*(ii-1)+1:ii*(nw+2));
        accelw(:,:,ii) = accelW(:,(nw+2)*(ii-1)+1:ii*(nw+2));
    end


    %save data every few timesteps
    if mod(i,N/framenumber)==0
        W(:,:,1+framenumber*i/N) = V0;
        d1wList(:,:,1+framenumber*i/N) = d1w;

        %estimate remaining time
        LoopTime1 = toc; LoopTime = (LoopTime1-LoopTime0) * framenumber/N * (N-i);
        fprintf('\n Estimated time remaining: %0.0f s\n', LoopTime)        
        LoopTime0 = LoopTime1;
        
    end
end
% % end main loop
% LoopTime = toc;
% LoopTime = LoopTime / N;

%save data
outputStruct = struct("W",W);
save(strcat(fileLabel,'.mat'), '-fromstruct', outputStruct);

%make video of simulation

%%%%%create AVI object
nFrames = size(W,3);
vidObj = VideoWriter(char(strcat(fileLabel)), 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);

colormap qualitative12

tubePlot = 0.6*tubeL;
diskRPlot = [diskR(1), diskR((1:nFrames-1)*N/framenumber)];
%%%%create movie
for j=1:1:nFrames
    clf;

    plotTangle(W(:,:,j),h, 'color',1:worms,'background','off');
    plotCylinder(diskRPlot(j)+h,tubePlot(j)+h, 'FaceAlpha',0.4);
    axis([(diskRPlot(j)+h)*[-1,1,-1,1],-tubePlot(j)-h,1*tubePlot(j)+h]);
    view(3); %caxis([0,11]);

    writeVideo(vidObj, getframe(gcf));
    drawnow;
    hold off;
    fprintf('\nCompiling video %0.0f%%\n', 100*j/nFrames);

end
close(vidObj);
hold off;





%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L0w, Ls0w, vLengths_fric] = initialLengths(W0)
    nw = size(W0,2)-2; worms = size(W0,1)/3;
  
    L0w = zeros(worms, nw+1);
    Ls0w = zeros(worms, nw+2);
    for i = 1:worms
        w = W0(3*(i-1)+1:3*i,:);
        L0w(i,:) = vecnorm(w(:,2:end) - w(:,1:end-1));
        Ls0w(i,:) = Vdom(w);
    end
    
    vLengths_fric = Ls0w'; 
    vLengths_fric(1,:) = 2*vLengths_fric(1,:); vLengths_fric(end,:) = 2*vLengths_fric(end,:);
    
end



function [velW, accelW] = velUpdate(worms, dt, FtotW, velW, accelW, V0comp, iMw, vLengths_fric, massList, L0w, eta_eff, zeta_eff, skindrag,rho )

    nw = size(L0w,2)-1;
    RvW = distsV(V0comp); %inverse squared distance
    iMcW = padarray(iMw, [1,1], 'post') + padarray(iMw,[1,1],'pre') - padarray(iMw,[1,1],'post').*padarray(iMw,[1,1],'pre');
    pDepthW = RvW.*iMcW; 
    pDepthW = (pDepthW .* vLengths_fric(:)) .* vLengths_fric(:);
    pDepthW = (pDepthW - diag(sum(pDepthW)));
    v0 = zeros(1,worms*(nw+2)); v1 = zeros(1,worms*(nw+2));
    for ii=1:worms
       v0((ii-1)*(nw+2)+1:ii*(nw+2)) = [1/L0w(ii,1), 1./L0w(ii,1:end-1)+1./L0w(ii,2:end), 1./L0w(ii,end)];
       v1((ii-1)*(nw+2)+1:ii*(nw+2)) = [1./L0w(ii,:), 0];
    end
    v1 = v1(1:end-1);
    Lapl0W = diag(v0) - diag(v1,1) - diag(v1,-1);
    LaplW = (eta_eff * Lapl0W - zeta_eff *pDepthW); 


    %linear operator A to get new velocity from old velocity, v_old
    A = (1+0.5*dt*skindrag/rho) * eye((nw+2)*worms) + 0.5*dt*(1./massList).*LaplW;
    v_old =  velW' + 0.5*dt*(FtotW'.*(1./massList') + accelW');

    
    %new velocity from solution of sparse linear system
    A = sparse(A);
    velW = lsqminnorm(A, v_old)';

    accelW = 1./massList .* (FtotW - skindrag/rho * massList .* velW - (LaplW *velW')' );


    % A = skindrag_eff * eye((nw+2)*worms) + LaplW;
    % A = sparse(A);
    % velW = lsqminnorm(A, FtotW')';

end






function Fcont = contactForce(w, diskR, tubeL, h, bdryForceStrength, fStrength)

    Fcont = zeros(size(w));
    Vradial = vecnorm(w(1:2,:));
    VLength = abs(w(3,:));
    dirns = w(:, Vradial>diskR-h);
    dirns = -dirns./vecnorm(dirns);

    bdryForce = 2*bdryForceStrength*dirns;
    Fcont(:, Vradial>diskR-h) = Fcont(:, Vradial>diskR-h) + bdryForce;

    dirns = w(3, VLength > tubeL-h);
    dirns = -[0*dirns; 0*dirns; sign(dirns)];

    bdryForce = 2*bdryForceStrength*dirns;
    Fcont(:, VLength>tubeL-h) = Fcont(:, VLength>tubeL-h) + bdryForce;

    ds = mean(vecnorm(w(:,2:end) - w(:,1:end-1)));
    Fcont = fStrength*ds*Fcont; %contact force with inverse friction parameter and length element

end




function B = pTransport(X,v0) %X is the set of positions, 3 by n+2, %v0 is initial vector
    n=size(X);
    n=n(2)-2;
    Y=zeros(3,n+1);
    Y(:,1)=v0; %initial frame
    for i=2:n+1
        ti=X(:,i)-X(:,i-1);
        tf=X(:,i+1)-X(:,i);
        R = rot(ti,tf);
        Y(:,i)=R*Y(:,i-1);
    end
    B=Y;
end


function R = rot(ti,tf)
        theta=atan2(norm(cross(ti,tf)),dot(ti,tf)); %angle to rotate ti into tf
        if norm(cross(ti,tf))==0
            m=[0,0,0];
        else
            m=cross(ti,tf)./norm(cross(ti,tf)); %the sign of m comes from here
        end
        m=theta*m;
        A=[0,-m(3),m(2);m(3),0,-m(1);-m(2),m(1),0]; %infinitesimal rotation
        R=expm(A);
end


function d1w = initialMaterialFrame(W0)
    nw = size(W0,2)-2; worms = size(W0,1)/3;
    d3w = W0(:,2) - W0(:,1);
    d1w = zeros(3*worms,nw+1);
    for i = 1:worms
        w = d3w(3*(i-1)+1:3*i);
        d1 = cross(cross(w, [0,0,1]'), w); 
        if norm(d1)==0
            d1 = cross(cross(w, [0,0,1.01]'), w);
        end
        d1=d1./norm(d1);
        d1w(3*(i-1)+1:3*i,:) = pTransport(W0(3*(i-1)+1:3*i,:), d1);
    end
end


function d1New = fupdate(X0,X1,d1)
    d30=X0(:,2:end)-X0(:,1:end-1); d31=X1(:,2:end)-X1(:,1:end-1);
    ph=atan2(sqrt(sum(cross(d30,d31).^2,1)),dot(d30,d31));
    N=cross(d30,d31);
    N=N./sqrt(dot(N,N));
    Y=cos(ph).*d1+((1-cos(ph)).*dot(N,d1)).*N + (sin(ph).*cross(N,d1));  
    Y(isnan(Y))=d1(isnan(Y));
    d1New=Y;
end



function R = dists(X) %distance between i'th and j'th links
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    L2 = diag(Xm'*Xm);
    R = -2*(Xm'*Xm) + (L2 + L2');
    R = R.^0.5; 
end


function R = distsV(X) %inverse squared distance between i'th and j'th points
    L2 = diag(X'*X);
    R = -2*(X'*X) + (L2 + L2');
    R = 1./(R + eye(size(R))) - eye(size(R));
end


function iM = interactMat(R, h) %Set R(i,j)=0 for |i-j|<=10
    n = size(R,1)-1;
    R([1:n+2:end, n+1+1:n+2:end-1]) = 2*h; 
    R([2:n+2:end-(n+1)]) = 2*h; 
    iM = R<2*h; %iM(i,j) = 1 iff link i and j are interacting
end




function pMat = pressureThreshold(R,iM,h,Kbulk,adhStrength)
    R(iM==0) = 0;
    pMat = Kbulk*( 1*(1 - R./(2*h)) + 100*(1-R./(2*h)).^3 + 100*(1-R./(2*h)).^5 + 100*(1-R./(2*h)).^7); %set cubic term to around 100
    pMat(iM==0)=0;
    pMat(pMat<0) = -adhStrength;
end




function Fi = iForce(X, pMat, R)
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    R(R==0) = 1;
    pMat = pMat./R;
    pMat = sum(pMat).*Xm - (Xm*pMat);
    pMat=[pMat(:,1), 0.5*(pMat(:,2:end)+pMat(:,1:end-1)), pMat(:,end)];
    Fi=pMat;
end



function Fs = sForce(X0,L0,gamma)
    d3 = X0(:,2:end)-X0(:,1:end-1);
    Ln = vecnorm(d3); d3 = d3./Ln;
    stretch = Ln./L0 - 1;
    F = gamma.*(stretch + 0*stretch.^3).*d3;
    Fs = [F, [0,0,0]'] + [[0,0,0]', -F];
end



function Ls = Vdom(X) %vector containing size of voronoi domain at each X
    n=size(X,2)-2;
    Y=(X-[X(:,1),X(:,1:n+1)]);
    Z=([X(:,2:n+2),X(:,n+2)]-X);
    Ls=0.5*(sqrt(sum(Y.^2))+sqrt(sum(Z.^2)));
end


function derivs = dL(X) %derivs(i,j,:) = dL(i)/dx(j) where L(i) is twice the voronoi domain of i'th vertex
    n=size(X);
    n=n(2)-2;
    Yp=zeros(n+2,n+2,3);
    ld=(X(:,1:n+1)-X(:,2:n+2))./(sqrt(sum((X(:,1:n+1)-X(:,2:n+2)).^2)));
    ud=-ld;
    d=[[0,0,0]',ud]-[ud,[0,0,0]'];
    Yp(:,:,1)=diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1);
    Yp(:,:,2)=diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1); 
    Yp(:,:,3)=diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1);
    derivs=Yp;
end


function ph1 = curvatureAngle1(X,d1)
    d3 = X(:,2:end)-X(:,1:end-1); d3 = d3./vecnorm(d3);
    d2=cross(d3, d1);
    ph1 = -0.5* dot( d2(:,2:end) + d2(:,1:end-1), d3(:,2:end) - d3(:,1:end-1)); 
    ph1=[0,ph1,0];
end

function ph2 = curvatureAngle2(X,d1)
    d3 = X(:,2:end)-X(:,1:end-1); d3 = d3./vecnorm(d3);
    ph2 = 0.5* dot( d1(:,2:end) + d1(:,1:end-1), d3(:,2:end) - d3(:,1:end-1)); 
    ph2=[0,ph2,0];
end


function dp = dphi1(X, d1) %dp(i,j,:) = dphi1(i)/dx(j), phi1 is kappa1 ds
    n=size(X,2)-2;
    Yp=zeros(n+2,n+2,3);
    d3=X(:,2:n+2)-X(:,1:n+1);
    l = vecnorm(d3); t = d3./vecnorm(d3);
    
    dm = cross(t,d1);
    
    ld = (1./l(1:end-1)) .* ( dm(:,2:end) + dot(t(:,2:end),t(:,1:end-1)).*dm(:,1:end-1) - t(:,1:end-1) .* dot(t(:,1:end-1),dm(:,2:end)) );   
    ud = (1./l(2:end)) .* ( dm(:,1:end-1) + dot(t(:,2:end),t(:,1:end-1)).*dm(:,2:end) - t(:,2:end) .* dot(t(:,2:end),dm(:,1:end-1)) );
    d=[[0,0,0]',-(ld+ud),[0,0,0]'];
    ld=[ld,[0,0,0]']; ud=[[0,0,0]',ud];

    Yp(:,:,1)=(diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1));
    Yp(:,:,2)=(diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1));
    Yp(:,:,3)=(diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1));
    Yp=real(Yp);
    Yp(isnan(Yp))=0; Yp(isinf(Yp))=0;
    dp = -0.5*Yp; % -1/2 for k1
end


function dp = dphi2(X, d1) %dp(i,j,:) = dphi2(i)/dx(j), phi2 is kappa2 ds
    n=size(X,2)-2;
    Yp=zeros(n+2,n+2,3);
    d3=X(:,2:n+2)-X(:,1:n+1);
    l = vecnorm(d3); t = d3./vecnorm(d3);
    
    dm = d1;
    
    ld = (1./l(1:end-1)) .* ( dm(:,2:end) + dot(t(:,2:end),t(:,1:end-1)).*dm(:,1:end-1) - t(:,1:end-1) .* dot(t(:,1:end-1),dm(:,2:end)) );   
    ud = (1./l(2:end)) .* ( dm(:,1:end-1) + dot(t(:,2:end),t(:,1:end-1)).*dm(:,2:end) - t(:,2:end) .* dot(t(:,2:end),dm(:,1:end-1)) );
    d=[[0,0,0]',-(ld+ud),[0,0,0]'];
    ld=[ld,[0,0,0]']; ud=[[0,0,0]',ud];

    Yp(:,:,1)=(diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1));
    Yp(:,:,2)=(diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1));
    Yp(:,:,3)=(diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1));
    Yp=real(Yp);
    Yp(isnan(Yp))=0; Yp(isinf(Yp))=0;
    dp = 0.5*Yp; % +1/2 for k2
end



function Fb=bForce_aniso1(X, d1, prefK1, alpha)
    n=size(X,2) - 2;
    ph=curvatureAngle1(X, d1); dph=dphi1(X, d1);
    L=2*Vdom(X); dl=dL(X);
    k1 = [0, zeros(1,n)+prefK1, 0]/2; %factor of 1/2 since dl, L are double Voronoi domain
    Yb=zeros(3,n+2);
    Yb(1,:)= -((2*alpha*ph)./L)*dph(:,:,1) + ((alpha*ph.^2)./L.^2)*dl(:,:,1);
    Yb(2,:)= -((2*alpha*ph)./L)*dph(:,:,2) + ((alpha*ph.^2)./L.^2)*dl(:,:,2);
    Yb(3,:)= -((2*alpha*ph)./L)*dph(:,:,3) + ((alpha*ph.^2)./L.^2)*dl(:,:,3);

    Yk=zeros(3,n+2);
    % Yk(1,:)= (2*alpha*prefK1)*sum(dph(:,:,1)) - (alpha*prefK1.^2)*sum(dl(:,:,1));
    % Yk(2,:)= (2*alpha*prefK1)*sum(dph(:,:,2)) - (alpha*prefK1.^2)*sum(dl(:,:,2));
    % Yk(3,:)= (2*alpha*prefK1)*sum(dph(:,:,3)) - (alpha*prefK1.^2)*sum(dl(:,:,3));
    Yk(1,:)= (2*alpha*k1)*dph(:,:,1) - (alpha*k1.^2)*dl(:,:,1);
    Yk(2,:)= (2*alpha*k1)*dph(:,:,2) - (alpha*k1.^2)*dl(:,:,2);
    Yk(3,:)= (2*alpha*k1)*dph(:,:,3) - (alpha*k1.^2)*dl(:,:,3);

    Fb=Yb+Yk;
end


function Fb=bForce_aniso2(X, d1, prefK2, alpha)
    n=size(X,2)-2;
    ph=curvatureAngle2(X, d1); dph=dphi2(X, d1);
    L=2*Vdom(X); dl=dL(X);
    k2 = [0, zeros(1,n)+prefK2, 0]/2; %factor of 1/2 since dl, L are double Voronoi domain
    Yb=zeros(3,n+2);
    Yb(1,:)= -((2*alpha*ph)./L)*dph(:,:,1) + ((alpha*ph.^2)./L.^2)*dl(:,:,1);
    Yb(2,:)= -((2*alpha*ph)./L)*dph(:,:,2) + ((alpha*ph.^2)./L.^2)*dl(:,:,2);
    Yb(3,:)= -((2*alpha*ph)./L)*dph(:,:,3) + ((alpha*ph.^2)./L.^2)*dl(:,:,3);

    Yk=zeros(3,n+2);
    % Yk(1,:)= (2*alpha*prefK2)*sum(dph(:,:,1)) - (alpha*prefK2.^2)*sum(dl(:,:,1));
    % Yk(2,:)= (2*alpha*prefK2)*sum(dph(:,:,2)) - (alpha*prefK2.^2)*sum(dl(:,:,2));
    % Yk(3,:)= (2*alpha*prefK2)*sum(dph(:,:,3)) - (alpha*prefK2.^2)*sum(dl(:,:,3));

    Yk(1,:)= (2*alpha*k2)*dph(:,:,1) - (alpha*k2.^2)*dl(:,:,1);
    Yk(2,:)= (2*alpha*k2)*dph(:,:,2) - (alpha*k2.^2)*dl(:,:,2);
    Yk(3,:)= (2*alpha*k2)*dph(:,:,3) - (alpha*k2.^2)*dl(:,:,3);

    Fb=Yb+Yk;
end



function [x,y,z,C]=tubeplot2(curve,r,n,ct,S)
% Usage: same as above but this gives colours the rod according to a 
% scalar field S along the curve.

  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
     end
  end
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end
  if nargin<4 || isempty(ct)
    ct=0.5*r;
  end


  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end

  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [~,idx]=min(abs(dv(:,1))); nvec(idx)=1;

  xyz=zeros(3,n+1,npoints+2);
  Col=zeros(3,n+1,npoints+2); 

  %precalculate cos and sing factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);

  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
    Col(:,:,k+1)=S(k);
  end
  %finally, cap the ends:
  xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  Ct=squeeze(Col(1,:,:));
  C=Ct; C(:,1) = C(:,2); C(:,end) = C(:,end-1);
  %... and plot:
  if nargout<3, surf(x,y,z,C); end
end





function plotTangle(W, h, varargin)

    nw = size(W,2); worms = size(W,1)/3;
    S = zeros(worms, nw+2);
    background = "on";
    lightSource = "on";
    wormAlpha = 1+zeros(1,worms);
    axShape = 1.2*max(abs(W(:)))*[-1,1,-1,1,-1,1];
    
    optInputs = varargin;
    for i = 1:length(optInputs)
        if string(optInputs{i}) == "color"
            S = optInputs{i+1}; 
        end
        if string(optInputs{i}) == "background"
            background = optInputs{i+1}; 
        end
        if string(optInputs{i}) == "lighting"
            lightSource = optInputs{i+1}; 
        end
        if string(optInputs{i}) == "tubeAlpha"
            wormAlpha = optInputs{i+1}; 
        end
        if string(optInputs{i}) == "axis"
            axShape = optInputs{i+1}; 
        end
    end
    
    if size(S,1)*size(S,2) == worms
        if size(S,2) == worms
            S = S';
        end
        S2 = zeros(worms, nw+2);
        S2 = S2 + S;
        S = S2;
    end
    
    for ii = 1:worms
        [x,y,z,C] = tubeplot2(W(3*(ii-1)+1:3*ii,:), h, 20, h/5, S(ii,:));
        Aplot=surf(x,y,z,C);
        set(Aplot,'FaceAlpha',wormAlpha(ii));
        shading interp;
        set(Aplot,'meshstyle','row');
        set(Aplot, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
        material shiny; 
        hold on;
    end

    %set frame
    axis(axShape);

    view(3);
    if lightSource == "on"
        camlight('left');
    end

    daspect([1,1,1]); 
    grid off; box off; axis off;
    set(gca, 'XTick',[]); set(gca, 'YTick',[]); set(gca,'ZTick',[]);
    
    if background == "on"
        box on; axis on;
        set(gca,'Color',[0.6 0.6 0.6]);
    end
    
    set(gcf, 'Color', [1,1,1]);
    set(gca,'LooseInset',get(gca,'TightInset'));
end





function plotCylinder(diskR, tubeL, varargin)

    faceAlpha = 0.5;
    optInputs = varargin;
    for i = 1:length(optInputs)
        if string(optInputs{i}) == "FaceAlpha"
            faceAlpha = optInputs{i+1}; 
        end
    end


    [X,Y,Z] = cylinder(diskR,50); Z = Z*2*tubeL-tubeL; C = zeros([size(Z),3])+0.6;
    cPlot = surf(X,Y,Z,C);
    shading interp;
    set(cPlot,'meshstyle','row');
    set(cPlot, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
    set(cPlot,'FaceAlpha',faceAlpha);

    T = linspace(0,2*pi,50);
    cTop = patch(diskR*cos(T), diskR*sin(T), 0*T+tubeL, squeeze(C(1,1,:))' );
    cBot = patch(diskR*cos(T), diskR*sin(T), 0*T-tubeL, squeeze(C(1,1,:))' );
    cTop.FaceAlpha = faceAlpha; cBot.FaceAlpha = faceAlpha;
end