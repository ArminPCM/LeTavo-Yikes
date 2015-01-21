%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ci_X_cip1, m_H_ci] = f_Compute_Velocities(c0_X_ci, m_R_c0)
%% Function that estimates the relative motion (velocities) between camera poses.
% Input:
% m_X 6xn matrix with n poses of the camera. The first three columns
% correspond to tx, ty, tz, and the last three to (RZ)roll, (RY)pitch and (RX)yaw.
% m_R_c0: Rotation from m to c0
% Output:
% ci_X_cip1: 6x(n-1) matrix with the  relative motions

if nargin < 2 || isempty(m_R_c0),
    m_R_c0 = rotox(-pi/2); 
end
c0_rpy_ci = c0_X_ci(4:6, :); % rpy rotations with respect to c0
c0_t_ci = c0_X_ci(1:3, :); % translations with respect to ci and matlab frame
numFrames = size(c0_rpy_ci, 2);

m_t_ci = m_R_c0*c0_t_ci; % translations with respect to ci and matlab frame
m_R_ci = cell(numFrames, 1);
m_H_ci = cell(numFrames, 1);
ci_X_cip1 =zeros(6, numFrames-1);

c0_R_ci = f_rpy2R(c0_rpy_ci([1 3 2], 1));
m_R_ci(1) = {m_R_c0*c0_R_ci};
m_H_ci(1) = {[[m_R_ci{1} m_t_ci(:, 1)]; 0 0 0 1]};
for i_frame=2:numFrames,
    c0_R_ci = f_rpy2R(c0_rpy_ci([1 3 2],  i_frame));
    m_R_ci(i_frame) = {m_R_c0*c0_R_ci};
    ci_X_cip1(1:3, i_frame-1) = m_R_ci{i_frame-1}' * (m_t_ci(:, i_frame) - m_t_ci(:, i_frame-1));
    ci_R_cip1 = m_R_ci{i_frame-1}' * m_R_ci{i_frame};
    [r, y, p] = f_R2rpy(ci_R_cip1);
    ci_X_cip1(4:6, i_frame-1) = mod([r, p, y]+pi, 2*pi)-pi;
    m_H_ci(i_frame) = {[[m_R_ci{i_frame} m_t_ci(:, i_frame)]; 0 0 0 1]};
end