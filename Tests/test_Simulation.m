%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Script for Simulation Case %%
% G. Puerto-Souza
% Astra Lab

% Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Comments and TODO's
% The training/Testing system is broken, fix it!
% add titles to plots and windows as well

close all
clear
clc
% run matlabpool or parpool to reduce processing time. Then uncomment
% parfor's and comment corresponding for's.
paralel_option = 0;
if paralel_option,
    parpool; % replace for matlabpool for matlab12 or older
end

%% Training Stage %%

%% 0.a Parameters Descriptors
addpath(genpath('..\ANN_Code'));
OF_Datasets = {'..\PreComputedData\Simulation\OF\Raw_OF_Blender_DHMA', '..\PreComputedData\Simulation\OF\Raw_OF_Blender_HMA', '..\PreComputedData\Simulation\OF\Raw_OF_Blender_Brox', '..\PreComputedData\Simulation\OF\Raw_OF_Blender_LK'};
OF_Datasets_Names = {'DHMA', 'HMA', 'Brox', 'LK'}; % same order than above
FLAG_ComputeDescriptors = 0;
SAVING_ComputedDescriptors = '..\PreComputedData\Simulation\Descriptors\Descriptors_DHBL'; % please specify route
% Note, IndSelectedFrames indicate all used frames. This means that if
% IndSelectedFrames has length n, then the optical fllow vectors are n-1
IndSelectedFrames =  [9:709]; %all frames -> load it from dataset!
%% 0.b Parameters Trajectory
INPUT_trajectory_filename = '..\Inputs\Simulation\Simulation_camera_pose_reversed_downscaled10.txt';
Visualize_ReadTrajectory = 0;
num_OF_Datasets_Names = length(OF_Datasets_Names);
%% I - Extracting Descriptors
if FLAG_ComputeDescriptors,
    %parameters
    Parameters_Descriptor.EnablingGrid = 1;
    Parameters_Descriptor.GridSize = [5, 5];
    Parameters_Descriptor.ImageSize = [640, 480];
    Parameters_Descriptor.GridLimits = [20, 600, 15, 450];
    Parameters_Descriptor.Units = [];% 'cartesian';%[];%'polar';
    Parameters_Descriptor.DescriptorFunction = [];%@(x) (mean(x,2));
    % DescriptorParameters.DescriptorFunctionName = 'mean'; % only if defining new descriptorfunction
    Parameters_Descriptor.MinThreshold = 0;
    Parameters_Descriptor.Normalize = 1;
    Parameters_Descriptor.SigmaNoise = 0;
    Parameters_Descriptor.EnablingLumen = 0; % need data about lumen... incomplete for now
    % DescriptorParameters.LumenAngles = -pi:pi/2:pi;
    % DescriptorParameters.LumenLevels = [50 50 50 100 100 Inf];
    for i_method=1:num_OF_Datasets_Names,
        %% Parameters
        Parameters_Input.dataset = OF_Datasets{i_method};
        Parameters_Input.strmethod = OF_Datasets_Names{i_method};
        %
        Descriptors.(OF_Datasets_Names{i_method})= f_Generate_Descriptors(Parameters_Input, IndSelectedFrames(1:end-1), Parameters_Descriptor);
        save(SAVING_ComputedDescriptors, 'Descriptors'); % conservative measure.
    end
else
    load(SAVING_ComputedDescriptors, 'Descriptors'); % loading saved descriptors
end

%%% GP -> TO DO: here we should define the training and testing sets!!

%% II - Loading Trajectory and simulating-noise/filtering

Parameters_Trajectory.Visualize = Visualize_ReadTrajectory;
Parameters_Trajectory.Vis_dx = 35;
wsize = 3;
[m_H_c, c0_X, ci_X_cip1, c0_X_int, m_H_int] = f_Read_Trajectory_from_Blender_file(INPUT_trajectory_filename, IndSelectedFrames, Parameters_Trajectory);
% IndSelectedFrames2Filter = [9-wsize:709+wsize];
IndSelectedFrames2Filter = [IndSelectedFrames(1)-wsize:IndSelectedFrames(end)+wsize];
[~, c0_X2filter, ~, ~, ~] = f_Read_Trajectory_from_Blender_file(INPUT_trajectory_filename, IndSelectedFrames2Filter, Parameters_Trajectory);
distances = sqrt(sum(ci_X_cip1(1:3, :).^2));
total_distance = sum(distances);
percentagefactor = (100/total_distance);
ColonScaling = 1.5/total_distance;
%% Add noise to the trajectory

noiset = 0.0008/ColonScaling;
noisea = deg2rad(1);
c0_X_noisy1 = [];c0_X_noisy2 = [];
numFrames_aux = size(c0_X2filter,2);
c0_X_noisy1(1:3, :) = c0_X2filter(1:3, :) - noiset*randn(3, numFrames_aux);
c0_X_noisy1(4:6, :) = c0_X2filter(4:6, :) - noisea*randn(3, numFrames_aux);
%
c0_X_noisy2(1:3, :) = c0_X2filter(1:3, :) - noiset*randn(3, numFrames_aux);
c0_X_noisy2(4:6, :) = c0_X2filter(4:6, :) - noisea*randn(3, numFrames_aux);


[ci_X_cip1_noisy1, m_H_ci_noisy1] = f_Compute_Velocities(c0_X_noisy1);
[ci_X_cip1_noisy2, m_H_ci_noisy2] = f_Compute_Velocities(c0_X_noisy2);

%%% DEBUG %%% uncomment to plot the trajectory in 3D
% figHandle = f_Plot_Trajectory(m_H_c, 1:41:698, [], 0.75, 0.5,  'b', 'g');
% f_Plot_Trajectory(m_H_ci_noisy1, 1:41:698, figHandle, 0.75, 0.5,  'r', 'n');
% f_Plot_Trajectory(m_H_ci_noisy2, 1:41:698, figHandle, 0.75, 0.5,  'r', 'n');
% plotting GT and noise

% gaussian filter
sigma = 5.0;
sizeg = 2*wsize;
x = linspace(-sizeg / 2, sizeg / 2, sizeg);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
%
c0_X_filtered1 = cell2mat(cellfun(@(x)(conv(x, gaussFilter, 'same')), num2cell(c0_X_noisy1, 2),'UniformOutput', false));
c0_X_filtered2 = cell2mat(cellfun(@(x)(conv(x, gaussFilter, 'same')), num2cell(c0_X_noisy2, 2),'UniformOutput', false));
c0_X_filtered1(:, [1:wsize end-wsize+1:end]) = [];
c0_X_filtered2(:, [1:wsize end-wsize+1:end]) = [];
c0_X_noisy1(:, [1:wsize end-wsize+1:end]) = [];
c0_X_noisy2(:, [1:wsize end-wsize+1:end]) = [];

%%% DEBUG %%% Uncomment to plot each velocities (useful to standarize domain of rpy)
% for i=1:6,
%     figure, plot(c0_X(i, :), 'k-');hold on;
%     plot(c0_X_noisy1(i, :), 'r-');
%     plot(c0_X_filtered1(i, :), 'b-');
% end
% for i=1:6,
%     figure, plot(c0_X(i, :), 'k-');hold on;
%     plot(c0_X_noisy2(i, :), 'r-');
%     plot(c0_X_filtered2(i, :), 'b-');
% end
%
[ci_X_cip1_filtered1, m_H_ci_filtered1] = f_Compute_Velocities(c0_X_filtered1);
%%% DEBUG %%% uncomment to plot the trajectory in 3D
% figHandle = f_Plot_Trajectory(m_H_c, 1:41:698, [], 0.75, 0.5,  'b', 'g');
% f_Plot_Trajectory(m_H_ci_noisy1, 1:41:698, figHandle, 0.75, 0.5,  'r', 'n');
% f_Plot_Trajectory(m_H_ci_filtered1, 1:41:698, figHandle, 0.75, 0.5,  'g', 'f');
% %
[ci_X_cip1_filtered2, m_H_ci_filtered2] = f_Compute_Velocities(c0_X_filtered2);
%%% DEBUG %%% uncomment to plot the trajectory in 3D
% figHandle = f_Plot_Trajectory(m_H_c, 1:41:698, [], 0.75, 0.5,  'b', 'g');
% f_Plot_Trajectory(m_H_ci_noisy2, 1:41:698, figHandle, 0.75, 0.5,  'r', 'n');
% f_Plot_Trajectory(m_H_ci_filtered2, 1:41:698, figHandle, 0.75, 0.5,  'g', 'f');

%% III - ANN Training %%
Parameters_ANNTraining.type_seeds = 'uniform'; %options: 'uniform', 'random' ... any other idea?
Parameters_ANNTraining.num_seeds = 100;
Parameters_ANNTraining.train_sz = 75;
Parameters_ANNTraining.val_sz = 15;
Parameters_ANNTraining.test_sz = 10;
Parameters_ANNTraining.seeds = f_Generate_Seeds(Parameters_ANNTraining.type_seeds, Parameters_ANNTraining.num_seeds);
% Descriptors2Use = []; %all descriptors, otherwise specify, e.g.:{'Grid_5x5_polar_mean', 'Grid_5x5_cartesian_mean', 'Grid_5x5_polar_median', 'Grid_5x5_cartesian_median'};
Descriptors2Use = {'Grid_5x5_polar_mean', 'Grid_5x5_cartesian_mean', 'Grid_5x5_polar_median', 'Grid_5x5_cartesian_median'};

%% Preparing Datasets
%training
ci_X_cip1_training = ci_X_cip1_filtered1;
Descriptors_training = Descriptors;
%testing
ci_X_cip1_testing = ci_X_cip1_filtered2;
Descriptors_testing = Descriptors;
% cells containing the ANN for Rot and Trans. We used cells in order to
% take advantage of the parallel toolbox
% Centering data
mean_ci_X_cip1_training = mean(ci_X_cip1_training, 2); % rowwise
ci_X_cip1_training = ci_X_cip1_training - mean_ci_X_cip1_training*ones(1, size( ci_X_cip1_training, 2));

ANNs_Translation = cell(num_OF_Datasets_Names, 1);
ANNs_Rotation = cell(num_OF_Datasets_Names, 1);
mean_Descriptors_training = cell(num_OF_Datasets_Names, 1);
% parfor inside!
for i_method=1:num_OF_Datasets_Names,
    [ANNs_Translation{i_method}, ANNs_Rotation{i_method}, mean_Descriptors_training{i_method}] = f_Train_ANN(ci_X_cip1_training, Descriptors_training.(OF_Datasets_Names{i_method}), Parameters_ANNTraining, Descriptors2Use, paralel_option);
end

%% IV - Predicting with the Trained ANN %%
% here!! extracting/loading descriptors for training
% for now assuming the same sets for training and testing.
c0_X_testing_hat = cell(num_OF_Datasets_Names, 1);
ci_X_cip1_testing_hat = cell(num_OF_Datasets_Names, 1);
m_H_testing_hat = cell(num_OF_Datasets_Names, 1);
ANN_Errors = cell(num_OF_Datasets_Names, 1);
% Predicting
for i_method=1:num_OF_Datasets_Names,
    [c0_X_testing_hat{i_method}, ci_X_cip1_testing_hat{i_method}, m_H_testing_hat{i_method}] = f_ANN_Prediction(Descriptors_testing.(OF_Datasets_Names{i_method}), mean_Descriptors_training{i_method},  mean_ci_X_cip1_training, ANNs_Translation{i_method}, ANNs_Rotation{i_method}, Parameters_ANNTraining, Descriptors2Use, paralel_option);
end
%% Plotting
% Computing best/median/wrost trajectories
controlPts = size(c0_X, 2);
for i_method=1:num_OF_Datasets_Names,
    ANN_Errors{i_method} = f_ANN_Plotting_Results(c0_X_testing_hat{i_method}, ci_X_cip1_testing_hat{i_method}, m_H_c, m_H_testing_hat{i_method}, c0_X, ci_X_cip1, controlPts, Parameters_ANNTraining, Descriptors2Use);
end

%% V - Final Plots %%
% Computing the Metrics
ind_NN2Plot = 1; % 1 plotting the best seeds, and 2 includes the median, and 3 includes the worst
E_CtrlPts =  cell(ind_NN2Plot, 1);
E_IntMax =  cell(ind_NN2Plot, 1);
E_VTran =  cell(ind_NN2Plot, 1);
E_VRot =  cell(ind_NN2Plot, 1);
numFramesVel = size(ci_X_cip1, 2);
for i_method=1:num_OF_Datasets_Names,
    for i_descriptors=1:length(Descriptors2Use),
        for i_seed = 1:ind_NN2Plot,
            ind2plot = ANN_Errors{i_method}.indices.(Descriptors2Use{i_descriptors})(i_seed);
            [E_CtrlPts{i_seed}, E_IntMax{i_seed}, E_VTran{i_seed}, E_VRot{i_seed}] = f_Compute_Metrics(...
                [], ci_X_cip1_testing_hat{i_method}.(Descriptors2Use{i_descriptors}){ind2plot}, m_H_int, ci_X_cip1, controlPts);
        end
        Metrics_ControlPts.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_CtrlPts;
        Metrics_IntTrajMax.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_IntMax;
        Metrics_Velocity_Translation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_VTran;
        Metrics_Velocity_Rotation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_VRot;
    end
end

% plotting - automatic saving was removed. TODO: add it again!
% 1 - Method/Descriptor Comparison
method_color = 'rbmk';
seed_symbol = 'v+^';
maxInt_symbol = 'hpo';
Descriptors2Use_plotnames = {'PM', 'CM', 'Pd', 'Cd'};
numDescrip = length(Descriptors2Use);
%% Control Pts
color_i = '  ';
figure('name','Control Points Error'); hold on;
title([' Control Points error for all methods']);
ylabel('Percentage error')
xlabel('Descriptor')
%creating labels
h = zeros(1, num_OF_Datasets_Names);
for i=1:num_OF_Datasets_Names,% descriptor
    h(i) = plot(0,0, method_color(i));
end
legend(OF_Datasets_Names(1:num_OF_Datasets_Names));
for i=1:num_OF_Datasets_Names,% descriptor
    set(h(i), 'Visible', 'off');
end

for i_method=1:num_OF_Datasets_Names, %method
    color_i(1) = method_color(i_method);
    for i_descriptors=1:numDescrip,
        x = numDescrip*(i_method-1) + i_descriptors;
        for i_seed = 1:ind_NN2Plot,
            % plotting all seeds specified (best/med/worst)
            y = percentagefactor*Metrics_ControlPts.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}){i_seed};
            color_i(2) = seed_symbol(i_seed);
            plot(x, y, color_i);
        end
    end
end
set(gca, 'Xtick', [1:numDescrip*num_OF_Datasets_Names])
set(gca, 'XtickL', repmat(Descriptors2Use_plotnames(1:numDescrip), 1,num_OF_Datasets_Names))
%% Max Int
color_i = '  ';
figure('name','Control Points Error'); hold on;
title([' Control Points error for all methods']);
ylabel('Percentage error')
xlabel('Descriptor')
%creating labels
h = zeros(1, num_OF_Datasets_Names);
for i=1:num_OF_Datasets_Names,% descriptor
    h(i) = plot(0,0, method_color(i));
end
legend(OF_Datasets_Names(1:num_OF_Datasets_Names));
for i=1:num_OF_Datasets_Names,% descriptor
    set(h(i), 'Visible', 'off');
end

for i_method=1:num_OF_Datasets_Names, %method
    color_i(1) = method_color(i_method);
    for i_descriptors=1:numDescrip,
        x = numDescrip*(i_method-1) + i_descriptors;
        for i_seed = 1:ind_NN2Plot,
            % plotting all seeds specified (best/med/worst)
            y = ColonScaling*Metrics_IntTrajMax.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}){i_seed};
            color_i(2) = maxInt_symbol(i_seed);
            plot(x, y, color_i);
        end
    end
end
set(gca, 'Xtick', [1:numDescrip*num_OF_Datasets_Names])
set(gca, 'XtickL', repmat(Descriptors2Use_plotnames(1:numDescrip), 1,num_OF_Datasets_Names))


% 2 - plotting the best translation/Rotation for each method
% since plotting several trajectories results very clutter, it is only done
% for the first one
for i_method=1:num_OF_Datasets_Names, %method
    figure('name',['Comparison of Velocity Error for ' OF_Datasets_Names{i_method}]);
    for i_descriptors=1:numDescrip,% descriptor
        subplot(2, numDescrip, i_descriptors);
        % plotting translation
        bestseed = 1; %median is 2 and worst is 3
        plot(ColonScaling*Metrics_Velocity_Translation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}){bestseed}, 'b-');
        xlabel([Descriptors2Use_plotnames{i_descriptors}, ' - time']);
        ylabel('$\| v-\widehat{v}\|$','interpreter','latex');
        ax = axis();
        axis([0 numFramesVel, ax(3:4)]);
        subplot(2, numDescrip, numDescrip+i_descriptors);
        %         plot(1-Metrics_Velocity_Rotation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}){bestseed}, 'r-');
        plot(Metrics_Velocity_Rotation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}){bestseed}, 'r-');
        xlabel([Descriptors2Use_plotnames{i_descriptors}, ' - time']);
        %         ylabel('$0.5(1-trace(R_{\omega}^T\widehat{R_{\omega}}))$','interpreter','latex');
        ylabel('degrees','interpreter','latex');
        ax = axis();
        axis([0 numFramesVel, ax(3:4)]);
    end
end
% plotting trajectories
for i_method=1:num_OF_Datasets_Names, %method
    figure('name',['Comparison of Integrated Trajectories for ' OF_Datasets_Names{i_method}]);
    for i_descriptors=1:numDescrip,% descriptor
        h = subplot(2, 2, i_descriptors);
        %h2 = f_Plot_Trajectory(m_H_int, 1:35:700, h, 0.75, 0.5,  'b', 'gt');
        %h2 = f_Plot_Trajectory(m_H_int, 1:35:700, h, 0.75, 0.5,  'b', []);
        h2 = f_Plot_Trajectory_LineOnly(m_H_int, 1:701, h, 'b-');
        indmin = ANN_Errors{i_method}.indices.(Descriptors2Use{i_descriptors})(1);
        %f_Plot_Trajectory(m_H_hat{i_method}.(Descriptors2Use{i_descriptors}){indmin}, 1:35:700, h, 0.75, 0.5,  'r', []);
        f_Plot_Trajectory_LineOnly(m_H_testing_hat{i_method}.(Descriptors2Use{i_descriptors}){indmin}, 1:701, h, 'r-');
        xlabel([Descriptors2Use_plotnames{i_descriptors}, ': b=GT, r=est']);
    end
end

%% Boxplots of seeds - ControlPts
figure;
for i_method=1:num_OF_Datasets_Names, %method
    x = zeros(Parameters_ANNTraining.num_seeds, numDescrip);
    for i_descriptors=1:numDescrip,
        x(:, i_descriptors) = ColonScaling*ANN_Errors{i_method}.ControlPts.(Descriptors2Use{i_descriptors});
    end
    subplot(2, 2, i_method);
    h = boxplot(x);
    set(gca, 'XTick', 1:numDescrip);
    set(gca, 'XTickLabel', Descriptors2Use_plotnames(1:numDescrip));
    xlabel(OF_Datasets_Names{i_method})
end

% all plotting modalities
figHandle = f_Plot_Trajectory_LineOnly_Animated(m_H_int, 1:701, [],'k');
figHandle = f_Plot_Trajectory_LineOnly(m_H_int, 1:701, [],'k');
% figHandle = f_Plot_Trajectory_LineOnly_Animated(m_H_testing_hat{3}.Grid_5x5_cartesian_median{1}, 1:701, figHandle,'r');
figHandle = f_Plot_Trajectory_LineOnly_Animated(m_H_testing_hat{2}.Grid_5x5_cartesian_median{2}, 1:701, figHandle,'g');
% figHandle = f_Plot_Trajectory_LineOnly_Animated(m_H_testing_hat{3}.Grid_5x5_cartesian_median{3}, 1:701, figHandle,'b');
% figHandle = f_Plot_Trajectory_LineOnly_Animated(m_H_testing_hat{2}.Grid_5x5_cartesian_median{4}, 1:701, figHandle,'m');
