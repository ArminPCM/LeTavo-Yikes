%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Error_ControlPts, Error_IntTrajMax, Error_Velocity_Translation, Error_Velocity_Rotation] = f_Compute_Metrics(c0_H_ci_hat, ci_X_cip1_hat, m_H_ci_GT, ci_X_cip1_GT, ControlPts, m_R_c0)
% plotting the performance of the pose prediction method
% output -> control points error, max integrated error, vel translation
% error, vel rotation error
%% Inputs
% c_H_ci_hat = estimated
% ControlPts = set of points in the integrated trajectory used to compute
% the error in the camera pose
% c_H_ci_GT =

if nargin < 7 || isempty(m_R_c0),
    m_R_c0 = rotox(-pi/2);
end

if ~isempty(ci_X_cip1_hat),
    [c0_X_int, m_H_ci_hat] = f_Integrated_Trajectory(ci_X_cip1_hat, m_R_c0);
else
    if ~isempty(c0_H_ci_hat),
        num = length(c0_H_ci_hat);
        m_H_c0 = eye(4);
        m_H_c0(1:3,1:3) = m_R_c0;
        for i=1:num,
            m_H_ci_hat(i) = {m_H_c0*c0_H_ci_hat{i}};
        end
        ci_X_cip1_hat = zeros(6, num-1);
        for i=1:num-1,
            ci_H_cip1_hat = c0_H_ci_hat{i}\c0_H_ci_hat{i+1};
            ci_X_cip1_hat(1:3, i) = ci_H_cip1_hat(1:3, 4);
            [r, y, p] = f_R2rpy(ci_H_cip1_hat(1:3, 1:3));
            ci_X_cip1_hat(4:6, i) = mod([r, p, y]+pi, 2*pi)-pi;
        end
    else
        error('velocities and poses cannot be empty at the same time')
    end
end


%checking error over control points
m_H_ci_hat_aux = [m_H_ci_hat{:}];
m_t_ci_hat = m_H_ci_hat_aux(1:3, 4:4:end);
m_H_ci_aux = [m_H_ci_GT{:}];
m_t_ci_GT = m_H_ci_aux(1:3, 4:4:end);
Diff_IntTrans = m_t_ci_GT - m_t_ci_hat;
Norm_ErrorTrans  = sqrt(sum(Diff_IntTrans.^2));
%% Translation over the Control Points
Error_ControlPts = sum(Norm_ErrorTrans(ControlPts));
%% now computing the maximum integrated trajectory error
Error_IntTrajMax = max(Norm_ErrorTrans);
%% Velocity Errors - Translation
Error_Velocity_Translation = sqrt(sum ((ci_X_cip1_hat(1:3, :) - ci_X_cip1_GT(1:3, :)).^2));
%% Velocity Errors - Rotation
% for i=1:size(ci_X_cip1_GT,2),
% 	ci_rpy_cip1 = ci_X_cip1_GT(4:6, i); % relative rpy rotation
%     ci_R_cip1 = f_rpy2R(ci_rpy_cip1([1 3 2]));
%     ci_rpy_cip1_hat = ci_X_cip1_hat(4:6, i); % relative rpy rotation
%     ci_R_cip1_hat = f_rpy2R(ci_rpy_cip1_hat([1 3 2]));
%     Error_Velocity_Rotation(i) = 0.5*trace(ci_R_cip1'*ci_R_cip1_hat) - 0.5;
% end
diff = mod(ci_X_cip1_hat(4:6, :) - ci_X_cip1_GT(4:6, :)+pi, 2*pi)-pi;
Error_Velocity_Rotation = rad2deg(sqrt(sum ((diff).^2)));

