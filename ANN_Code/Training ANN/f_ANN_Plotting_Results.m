%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ANN_Errors = f_ANN_Plotting_Results(c0_X_hat, ci_X_cip1_hat, m_H,  m_H_hat, GT_c0_X, GT_ci_X_cip1, controlPts, Parameters, Descriptors2Use)
%%
if nargin<9 || isempty(Descriptors2Use),
    Descriptors2Use = fieldnames(c0_X_hat);
end
display('Evaluating Performance!')
% metrics are the ending point distance, maximum velocity  errors, sum
% absolute erorr, absolute error
GT_controlPts = GT_c0_X(1:3, controlPts); %or over control points!
GT_c0_X_1_3 = GT_c0_X(1:3,:);
numVels = size(GT_ci_X_cip1, 2);
numctrlpts = size(GT_controlPts, 2);
Errors_normcontrolPts_i = zeros(Parameters.num_seeds, 1);
Errors_controlPts_i = zeros(numctrlpts, Parameters.num_seeds);
% Errors_MaxIntTraj_i =  zeros(6, Parameters.num_seeds);
for i_descriptor = 1:length(Descriptors2Use),
    display(Descriptors2Use{i_descriptor});
    ci_X_cip1_hat_i = ci_X_cip1_hat.(Descriptors2Use{i_descriptor});
    c0_X_hat_i = c0_X_hat.(Descriptors2Use{i_descriptor});
    for i_seed = 1:Parameters.num_seeds,
        c0_X_hat_i_ctrlPts = c0_X_hat_i{i_seed}(1:3, controlPts);
        %         c0_X_hat_i_1_3 = c0_X_hat_i{i_seed}(1:3, :);
        %% Integrated Trajectories
        % wrt the control points
        Errors_controlPts_i(:, i_seed) = sqrt(sum((GT_controlPts - c0_X_hat_i_ctrlPts).^2));
        Errors_normcontrolPts_i(i_seed) = sum(Errors_controlPts_i(:, i_seed));
        % maximum integrated error
        %         Error_IntegratedTrajectories_i = sqrt(sum((GT_c0_X_1_3 - c0_X_hat_i_1_3).^2));
        %         [maxValue, maxind] = max(Error_IntegratedTrajectories_i);
        %         Errors_MaxIntTraj_i(i_seed) = maxValue;
        % Velocities
        %         Diff_RelativeMotion = GT_ci_X_cip1 - ci_X_cip1_hat_i{i_seed};
        %         Errors_Velocity_i(:, i_seed) = maxValues;
        %         residualsT = sqrt(sum(Diff_RelativeMotion(1:3, :).^2));%norm of the errors
        %         Errors_meanResiduals_i(i_seed) = mean(residuals);
        %         Errors_TranslationalVelocity_aux(i_seed) = max(residuals);
    end
    % selecting the smallest control points error
    [minValue, indmin] = min(Errors_normcontrolPts_i);
    [sortedValue, indsort] = sort(Errors_normcontrolPts_i);
    indmed = indsort(floor(Parameters.num_seeds/2));
    [maxValue, indmax] = max(Errors_normcontrolPts_i);
    % computing velocity errors for the best seed
    indices = [indmin, indmed, indmax];
    values = [minValue, sortedValue(indmed), maxValue];
%     Errors_Velocity_i = cell(3, 1);
%     residualsT = zeros(3, numVels);
%     residualsR = zeros(3, numVels);
%     Errors_MaxIntTraj_i = zeros(3, 1);
%     for i_ind=1:3,
%         ind_curr = indices(i_ind);
%         Errors_Velocity_i(i_ind) = {GT_ci_X_cip1 - ci_X_cip1_hat_i{ind_curr}};
%         residualsT(i_ind, :) = sqrt(sum(Errors_Velocity_i{i_ind}(1:3, :).^2));%norm of the errors
%         for i_frames=1:numVels,
%             ci_rpy_cip1 = GT_ci_X_cip1(4:6, i_frames); % relative rpy rotation
%             ci_R_cip1 = f_rpy2R(ci_rpy_cip1([1 3 2]));
%             ci_rpy_cip1_hat = ci_X_cip1_hat_i{ind_curr}(4:6, i_frames); % relative rpy rotation
%             ci_R_cip1_hat = f_rpy2R(ci_rpy_cip1_hat([1 3 2]));
%             residualsR(i_ind, i_frames) = 0.5*trace(ci_R_cip1'*ci_R_cip1_hat) - 0.5;
%         end
%         c0_X_hat_i_1_3 = c0_X_hat_i{indices(i_ind)}(1:3, :);
%         Error_IntegratedTrajectories_i = sqrt(sum((GT_c0_X_1_3 - c0_X_hat_i_1_3).^2));
%         [maxValue, dummy] = max(Error_IntegratedTrajectories_i);
%         Errors_MaxIntTraj_i(i_ind) = maxValue;
%         if maxValue < minValue,
%             keyboard;
%         end
%     end
    Errors_controlPts_vectors.(Descriptors2Use{i_descriptor}) = Errors_controlPts_i;
    Errors_controlPts.(Descriptors2Use{i_descriptor}) = Errors_normcontrolPts_i;
%     Errors_MaxIntTraj.(Descriptors2Use{i_descriptor}) = Errors_MaxIntTraj_i;
%     Errors_Velocity.(Descriptors2Use{i_descriptor}) = Errors_Velocity_i;
%     Errors_VelocityTranslation.(Descriptors2Use{i_descriptor}) = residualsT;
%     Errors_VelocityRotation.(Descriptors2Use{i_descriptor}) = residualsR;
    Errors_indices.(Descriptors2Use{i_descriptor}) = indices;
    Errors_indvalues.(Descriptors2Use{i_descriptor}) = values;
%     plotting trajectories
%         figHandle = f_Plot_Trajectory(m_H, 1:35:700, [], 0.75, 0.5,  'b', 'gt');
%         f_Plot_Trajectory(m_H_hat.(Descriptors2Use{i_descriptor}){indmin}, 1:35:700, figHandle, 0.75, 0.5,  'r', ' ');
%         title(Descriptors2Use{i_descriptor});
%         for i_vel=1:6,
%             figure;
%             plot(GT_ci_X_cip1(i_vel, :), 'b-'); hold on
%             plot(ci_X_cip1_hat_i{indmin}(i_vel, :), 'r:');
%         end
%     
end
ANN_Errors.ControlPts = Errors_controlPts;
ANN_Errors.ControlPts_vectors = Errors_controlPts_vectors;
% ANN_Errors.MaxIntTraj = Errors_MaxIntTraj;
ANN_Errors.indices = Errors_indices;
ANN_Errors.indvalues = Errors_indvalues;
% ANN_Errors.Velocity = Errors_Velocity;
% ANN_Errors.VelocityTranslation = Errors_VelocityTranslation;
% ANN_Errors.VelocityRotation = Errors_VelocityRotation;
% ANN_Errors = {Errors_controlPts, Errors_controlPts_vectors, Errors_MaxIntTraj, Errors_indmin, Errors_Velocity, Errors_VelocityTranslation, Errors_VelocityRotation};