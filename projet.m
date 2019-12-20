Yini = single(imread('Mars_dunes.jpg')); 
ltot = size(Yini, 1);
ctot = size(Yini, 2);
trois = size(Yini,3); 
image(uint8(Yini));
title('image initiale');
axis equal

nl = 500;
global l;
l = floor(ltot/nl);
nc = 1000;
global c;
c = floor(ctot/nc);

disp("Avant Extraction des blocs");
A = zeros(nl * nc, l*c, 3); 
for i = 1:nl
    %disp("Extraction i des blocs");
    for j = 1:nc
        Y = Yini((i-1)*l+1: i*l, (j-1)*c+1:j*c,:);
        A((i-1)*nc + j,:, :) = reshape(Y, l*c, 3);
    end
end

disp("Apres Extraction des blocs");


Yfin = zeros(ltot, ctot, 3);

for i = 1:nl
    for j = 1:nc
        Yfin((i-1)*l+1: i*l, (j-1)*c+1:j*c,:) = reshape(A((i-1)*nc + j,:,:), l, c, 3);
    end
end

figure(2)
image(uint8(Yfin));
title('reconstruction test');
axis equal

Acompresse = zeros(nl * nc, l * c, 3);
p = 2;  %doit etre inferieur a 3
disp("Avant compression des blocs");
Iptot = 0;
for i = 1: nc * nl
    %disp("Avant codage des blocs");
    Atemp = A(i,:,:);
    Atemp = reshape(Atemp, l*c, 3);
    [P, E, Ip] = codeur_ACP(Atemp,p);
    %disp("Avant decodage des blocs");
    Atemp = decodeur_ACP(P,E);
    Atemp = reshape(Atemp,l*c, 3);
    Acompresse(i,:,:) = Atemp;
    Iptot = Iptot + Ip;
end

disp("A compresse contient ");
disp(numel(Acompresse));

Iptot = Iptot / (nc* nl);
disp("Apres compression des blocs");
disp("Iptot vaut")
disp(Iptot);

disp("Avant affichage");
for i = 1:nl
%     disp("dans affichage");
    for j = 1:nc
        Ycompresse((i-1)*l+1: i*l, (j-1)*c+1:j*c,:) = reshape(Acompresse((i-1)*nc + j,:,:), l, c, 3);
    end
end
disp("Apres affichage");

figure(3)
image(uint8(Ycompresse));
title('reconstruction compressee');
axis equal


% S1 = std(X(:,1))
% S2 = std(X(:,2))
% S3 = std(X(:,3))
% V = cov(X)
% [E,D] = eig(V)
% 
% P = X*E(:,2:3)
% Einv = E'
% X2 = P*Einv(2:3,:)
% Yfin = reshape(X2, ltot, ctot, 3);
% figure(10)
% image(uint8(Yfin));
% 
% 
% S4 = std(X2(:,1))
% S5 = std(X2(:,2))
% S6 = std(X2(:,3))
% 
% S = S1+S2+S3
% S10 = S4+S5+S6

function [P, E, Ip] = codeur_ACP(X,p)
    Xcov = cov(X);
    [E,D] = eig(Xcov);
    E = E(:,3-p+1: 3);
    P = X*E;
    Iptemp(1) = std(X(:,1));
    Iptemp(2) = std(X(:,2));
    Iptemp(3) = std(X(:,3));
    Ip = 0;
    for i = 3-p+1: 3
        Ip = Ip +  Iptemp(i);
    end
end

function X = decodeur_ACP(P,E)
    X = P * E';
end

