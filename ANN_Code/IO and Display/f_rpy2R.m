%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G-L. Mariottini
%  Astra Lab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R=f_rpy2R(V) 
%z axis is the front of the camera
roll=V(1);
pitch=V(2);
yaw=V(3);
phi=roll;
theta=pitch;
psi=yaw;

c_phi=cos(phi);
c_theta=cos(theta);
c_psi=cos(psi);
s_phi=sin(phi);
s_theta=sin(theta);
s_psi=sin(psi);


R = [ c_phi*c_theta, c_phi*s_theta*s_psi-s_phi*c_psi, c_phi*s_theta*c_psi+s_phi*s_psi;
      s_phi*c_theta, s_phi*s_theta*s_psi+c_phi*c_psi, s_phi*s_theta*c_psi-c_phi*s_psi;
           -s_theta,   c_theta*s_psi,                c_theta*c_psi];
% roll=V(1);
% pitch=V(2);
% yaw=V(3);
% gamma=roll;
% beta=pitch;
% alpha=yaw;
% 
% c_gamma=cos(gamma);
% c_beta=cos(beta);
% c_alpha=cos(alpha);
% s_gamma=sin(gamma);
% s_beta=sin(beta);
% s_alpha=sin(alpha);
% 
% 
% R = [ c_alpha*c_beta, c_alpha*s_beta*s_gamma-s_alpha*c_gamma, c_alpha*s_beta*c_gamma+s_gamma*s_alpha;
%       s_alpha*c_beta, s_alpha*s_beta*s_gamma+c_alpha*c_gamma, s_alpha*s_beta*c_gamma-s_gamma*c_alpha;
%            -s_beta,   c_beta*s_gamma,                c_beta*c_gamma];