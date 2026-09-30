%Code to generate initial conditions for swarm linking

clear;


%parameters of IC
worms = 20; %target number of filaments
h = 0.3; %radius of filaments
bradius = 2.8; %filaments sampled from ball with this radius
nSteps = 6; %sets the number of discrete points on each filament
ds = 2*h; %discretization length of filaments


%Generate new initial condition (stored in W)
W = zeros(3*worms, 2*nSteps+1);

%sample initial filament
x = normrnd(0,1,3,1); x = x./vecnorm(x);
x = x*bradius*rand;
n = normrnd(0,1,3,1); n = n./vecnorm(n);
W(1:3,:) = x + ds*n.*[-nSteps:1:nSteps];

ii = 1; %counter for current number of filaments
while ii < worms

    %generate candidate next filament
    x = normrnd(0,1,3,1); x = x./vecnorm(x);
    x = x*bradius*rand;
    n = normrnd(0,1,3,1); n = n./vecnorm(n);
    w = x + ds*n.*[-nSteps:1:nSteps];
    
    Wtest = [W(1:3*ii,:) ; w];

    B  = pagetranspose(reshape(Wtest',[size(Wtest,2),3,size(Wtest,1)/3]));
    B = reshape(B,3,[]);

    R = realDists(B);
    R = R + kron( eye(ii+1), ones(2*nSteps+1, 2*nSteps+1) );

    %accept candidation filament if it doesn't intersect any current filament
    if min(R(:)) >= 1.25*2*h
        W(3*ii+1:3*ii+3,:) = w;
        ii = ii+1;
        fprintf('\n Filaments generated: %0.0f out of %0.0f \n', ii, worms);
    end

end


%save new initial condition
save(strcat('initial_conditions/IC_new.mat'), 'W');

%plot new initial condition
clf;
plotTangle(W,h,'color',1:worms); axis auto;




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


function R = realDists(X) %squared distance between i'th and j'th points
    L2 = diag(X'*X);
    R = -2*(X'*X) + (L2 + L2');
    R = R.^0.5;
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



