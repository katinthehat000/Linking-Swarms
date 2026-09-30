%Code to generate contact link for a simulation of a linking swarm
%output is a graph showing how contact link varies with time


clear;


fileLabel = 'results/output_sim.mat';

h=0.3; %radius
load(fileLabel, 'W');
cLink = zeros(1,size(W,3));
for ii = 1:size(W,3)
    cLink(ii) = cLinkPerWorm(W(:,:,ii), 1.25*h);
end

clf; hold on;
plot(0:size(W,3)-1, cLink, 'LineWidth',2);

axis square; box on;
xlabel('Time (T_0)'); ylabel('Total contact link per curve');
set(gca, 'FontSize', 16)











function cLink = cLinkPerWorm(W, h)
    worms = size(W,1)/3;
    %geometric adjacency matrix
    [adjMat, ~] = geoAdjMat(W, h, worms); %distMat <1 implies contact

    LkMat0 = LinkMatrix(W, worms);

    %contact link matrix
    cLinkMat = LkMat0.*(adjMat~=0);
    cLinkMat(1:worms+1:end) = 0;

    cLink = sum(abs(cLinkMat(:)))/worms;       
end


function [adjMat, distMat] = geoAdjMat(W, h, worms)

    nw = size(W,2)-2;
    V0comp = zeros(3, worms*(nw+2));
    for ii = 1:worms
        V0comp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = W(3*(ii-1)+1:3*ii,:);
    end
    Rw = realdistsV(V0comp);

    adjMat = zeros(worms,worms);
    distMat = zeros(worms,worms);
    for i = 1:worms
        for j = 1:worms      
            Rwsub = Rw( (i-1)*(nw+2)+1:i*(nw+2), (j-1)*(nw+2)+1:j*(nw+2) );
            iMsub = Rwsub < 2*h;
            adjMat(i,j) = sign(nnz(iMsub));
            distMat(i,j) = min(Rwsub(:));
            if i==j
                distMat(i,j)=0;
                adjMat(i,j)=0;
            end
        end
    end

    distMat = distMat./(2*h); %so distMat < 1 -> contact

end




function R = realdistsV(X) %distance between i'th and j'th points
    L2 = diag(X'*X);
    R = -2*(X'*X) + (L2 + L2');
    R = R.^0.5;
end




function Omega = Link(X1,X2)
    Omega = 0;
    for i = 1:size(X1,2)-1
        for j = 1:size(X2,2)-1
            n0 = X1(:,i) - X2(:,j); nt = X1(:,i) - X2(:,j+1);
            ns = X1(:,i+1) - X2(:,j); nst = X1(:,i+1) - X2(:,j+1);

            E0_s = cross(n0,ns); E0_s = E0_s./vecnorm(E0_s);
            Es_st = cross(ns,nst); Es_st = Es_st./vecnorm(Es_st);
            Est_t = cross(nst,nt); Est_t = Est_t./vecnorm(Est_t);
            Et_0 = cross(nt,n0); Et_0 = Et_0./vecnorm(Et_0);

            A = asin(dot(Et_0, E0_s)) + asin(dot(E0_s, Es_st)) + asin(dot(Es_st, Est_t)) + asin(dot(Est_t, Et_0));

            Omega = Omega + A/(4*pi) * sign( dot(n0, cross(X1(:,i+1)-X1(:,i), X2(:,j+1)-X2(:,j))) );

        end
    end
end


function LkMat = LinkMatrix(W0, worms)
    LkMat = zeros(worms,worms);
    for i = 2:worms
        for j = 1:i-1
            LkMat(i,j) = Link(W0(3*i-2:3*i,:), W0(3*j-2:3*j,:));
        end
    end
    LkMat = LkMat + LkMat';
end
