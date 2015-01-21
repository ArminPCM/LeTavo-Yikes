%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [c0_X, m_H_ci] = f_Integrated_Trajectory(ci_X_cip1, m_R_c0)
%% Function that estimates the integrated trajectory from relative motions.
% Input:
% ci_X_cip1 6xn matrix with n relative poses of the camera. The first three columns
% correspond to tx, ty, tz, and the last three to (RZ)roll, (RY)pitch and (RX)yaw.
% m_R_c0: Rotation from m to c0
% Output:
% M_Traj_ci: 6x(n+1) matrix with the camera poses in the c0 reference frame


if nargin < 2 || isempty(m_R_c0),
    m_R_c0 = rotox(-pi/2); 
end

f_assemble_H_matrix = @(R, t)([[R t]; 0 0 0 1]); % function that constructs the projection matrix H.

numFrames = size(ci_X_cip1, 2);
m_R_ci = cell(numFrames+1, 1); % matrix storing rotations from matlab to the camera i
m_t_ci = zeros(3, numFrames+1); % matrix storing translations from matlab to the camera i
m_H_ci = cell(numFrames+1, 1); % projection matrix H for each camera pose in the matlab reference frame
c0_angles_ci = zeros(3, numFrames+1); % matrix with the rpy for ci at reference frame c0
% Initial point is set as [0;0;0];
m_R_ci(1) = {m_R_c0};
m_H_ci(1) = {f_assemble_H_matrix(m_R_c0, m_t_ci(:, 1))};
for i_frames=1:numFrames, 
    ci_rpy_cip1 = ci_X_cip1(4:6, i_frames); % relative rpy rotation 
    ci_R_cip1 = f_rpy2R(ci_rpy_cip1([1 3 2]));
	m_R_cip1 = m_R_ci{i_frames}*ci_R_cip1; 
    ci_t_cip1 = ci_X_cip1(1:3, i_frames); % relative translation 
    m_t_ci(:, i_frames+1) = m_t_ci(:, i_frames) + m_R_ci{i_frames}*ci_t_cip1;
    [r, y, p] = f_R2rpy(m_R_c0'*m_R_cip1);
    c0_angles_ci(:, i_frames+1) = mod([r, p, y], 2*pi);
    m_R_ci(i_frames+1) = {m_R_cip1};
    m_H_ci(i_frames+1) = {f_assemble_H_matrix(m_R_ci{i_frames+1}, m_t_ci(:, i_frames+1))};
end
% packing integrated trajectory
c0_X = [m_R_c0'*m_t_ci; c0_angles_ci];
