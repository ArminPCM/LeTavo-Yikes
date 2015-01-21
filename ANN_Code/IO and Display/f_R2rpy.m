%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G-L. Mariottini
%  Astra Lab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,p,y]=f_R2rpy(R,sol)
if nargin==1,
    sol=0;
end

if sol==0,
    phi = f_atan(R(2,1),R(1,1));
    theta = f_atan(-R(3,1), sqrt(R(3,2)^2+R(3,3)^2));
    psi = f_atan(R(3,2), R(3,3));
else
    phi = f_atan(-R(2,1),-R(1,1));
    theta = f_atan(-R(3,1), -sqrt(R(3,2)^2+R(3,3)^2));
    psi = f_atan(-R(3,2), -R(3,3));
end

r=phi;
p=theta;
y=psi;