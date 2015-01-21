%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Script for Simulation Case %%
% G. Puerto-Souza
% Astra Lab

% Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Animal Experiments -> To add remaining functions and check code
% compatibility
close all
clear
clc
%matlabpool or parpool
Parameters_ANNTraining.num_seeds = 100;
%% Descriptors Extracted from standalone script
SAVING_ComputedDescriptors = {'Descriptors_Piggy1_DHBL', 'Descriptors_Piggy2_DHBL'};
OF_Datasets_Names = {'DHMA', 'HMA', 'Brox', 'LK'};
%% Test 1: Using same dataset, the first one
load(SAVING_ComputedDescriptors{1},'IndSelectedFrames');
ind_testing = [2277:2360];%[2112:2156]; % wrt image number! Opical flow from i->i+1
[ind_training, ind_Descriptors_training] = setdiff(IndSelectedFrames, ind_testing);
[ind_testing, ~, ind_Descriptors_testing] = intersect(IndSelectedFrames, ind_testing);
%% Parameters Trajectory
INPUT_trajectory_filename = {'./Piggy/NDI_Data1.txt', './Piggy/NDI_Data2.txt'};
Visualize_ReadTrajectory = 1;
num_OF_Datasets_Names = length(OF_Datasets_Names);
% Descriptors2Use = {'Grid_5x5_polar_mean', 'Grid_5x5_cartesian_mean', 'Grid_5x5_polar_median', 'Grid_5x5_cartesian_median'};
Descriptors2Use = {'Grid_5x5_cartesian_mean', 'Grid_5x5_cartesian_median'};
num_Descriptors2Use = length(Descriptors2Use);
%% Loading Descriptors!
load(SAVING_ComputedDescriptors{1}, 'Descriptors'); % loading saved descriptors
for i_method=1:num_OF_Datasets_Names,
    for i_descriptor = 1:num_Descriptors2Use,
        Descriptors_training.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptor}) = ...
            Descriptors.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptor})(:, ind_Descriptors_training);
        Descriptors_testing.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptor}) = ...
            Descriptors.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptor})(:, ind_Descriptors_testing);
    end
end

%% Loading Trajectory
% reading the file
A_trackerData = dlmread(INPUT_trajectory_filename{1});
M_H_A = f_aurora2Matlab(A_trackerData);

% f_Plot_Trajectory(M_H_A, 1550:25:2413, [],10,5, 'r', ''); 
% f_Plot_Trajectory_LineOnly_Animated(M_H_A, 1550:2413, [],'r'); 
% computing velocities
ColonScaling = 1/1000;%scale is mm
ci_X_cip1 = f_Compute_Velocities_from_H(M_H_A);
% cropping good indices
% units are meters scaling to cm
ci_X_cip1(1:3, :) = ColonScaling*ci_X_cip1(1:3, :);
ci_X_cip1_rawTrain = ci_X_cip1(:, ind_training);
ci_X_cip1_rawTest = ci_X_cip1(:, ind_testing);
% Integrating them
[c0_X_rawTrain,  m_H_rawTrain] = f_Integrated_Trajectory(ci_X_cip1_rawTrain); % integrating trajectories
c0_X_rawTrain(4, 6:13) = c0_X_rawTrain(4, 6:13)-2*pi;
c0_X_rawTrain(5, :) = mod(c0_X_rawTrain(5, :)+pi, 2*pi)-pi;
c0_X_rawTrain(6, :) = mod(c0_X_rawTrain(6, :)+pi/2, 2*pi)-pi/2;
% f_Plot_Trajectory(m_H_rawTrain, 1:10:length(m_H_rawTrain), [],10,5, 'r', ''); 
[c0_X_rawTest,  m_H_rawTest] = f_Integrated_Trajectory(ci_X_cip1_rawTest); % integrating trajectories
c0_X_rawTest(4, 2:4) = c0_X_rawTest(4, 2:4)-2*pi;
c0_X_rawTest(5, :) = mod(c0_X_rawTest(5, :)+pi, 2*pi)-pi;
c0_X_rawTest(6, :) = mod(c0_X_rawTest(6, :)+pi, 2*pi)-pi;


% f_Plot_Trajectory(m_H_rawTest, 1:10:length(m_H_rawTest), [],10,5, 'r', ''); 
% Add noise to the trajectory
% gaussian filter parameters
window_size = 2;
sigma = 0.500;
sizeg = 2*window_size;
x = linspace(-sizeg / 2, sizeg / 2, sizeg);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
% adding window to the original data to filter
% % mirror type
% c0_X_raw_4FilteringTrain = c0_X_rawTrain(:, [window_size:-1:1 1:end  end-window_size+1:end]);
% c0_X_raw_4FilteringTe4st = c0_X_rawTest(:, [window_size:-1:1 1:end  end-window_size+1:end]);
% same type
c0_X_raw_4FilteringTrain = c0_X_rawTrain(:, [ones(1, window_size) 1:end  end*ones(1, window_size)]);
c0_X_raw_4FilteringTest = c0_X_rawTest(:, [ones(1, window_size) 1:end  end*ones(1, window_size)]);

c0_X_filteredTrain = cell2mat(cellfun(@(x)(conv(x, gaussFilter, 'same')), num2cell(c0_X_raw_4FilteringTrain, 2),'UniformOutput', false));
c0_X_filteredTest = cell2mat(cellfun(@(x)(conv(x, gaussFilter, 'same')), num2cell(c0_X_raw_4FilteringTest, 2),'UniformOutput', false));
c0_X_filteredTrain(:, [1:window_size end-window_size+1:end]) = [];
c0_X_filteredTest(:, [1:window_size end-window_size+1:end]) = [];
c0_X_filteredTrain(:, 1) = zeros(6,1);
c0_X_filteredTest(:, 1) = zeros(6,1);

% for i=1:6,
%     figure, plot(c0_X_rawTrain(i, :), 'k-');hold on;
%     plot(c0_X_filteredTrain(i, :), 'b-');
% end
% 
% for i=1:6,
%     figure, plot(c0_X_rawTest(i, :), 'k-');hold on;
%     plot(c0_X_filteredTest(i, :), 'b-');
% end

[ci_X_cip1_filteredTrain, m_H_ci_filtetredTrain] = f_Compute_Velocities(c0_X_filteredTrain);
% figHandle = f_Plot_Trajectory(m_H_c, 1:41:698, [], 0.75, 0.5,  'b', 'g');
% f_Plot_Trajectory(m_H_ci_noisy1, 1:41:698, figHandle, 0.75, 0.5,  'r', 'n');
% f_Plot_Trajectory(m_H_ci_filtered1, 1:41:698, figHandle, 0.75, 0.5,  'g', 'f');
% %
[ci_X_cip1_filteredTest, m_H_ci_filteredTest] = f_Compute_Velocities(c0_X_filteredTest);
% figHandle = f_Plot_Trajectory(m_H_c, 1:41:698, [], 0.75, 0.5,  'b', 'g');
% f_Plot_Trajectory(m_H_ci_noisy2, 1:41:698, figHandle, 0.75, 0.5,  'r', 'n');
% f_Plot_Trajectory(m_H_ci_filtered2, 1:41:698, figHandle, 0.75, 0.5,  'g', 'f');
% 
% ci_X_cip1_filteredTrain = ci_X_cip1_rawTrain;
ci_X_cip1_filteredTest = ci_X_cip1_rawTest;
m_H_ci_filteredTest = m_H_rawTest;
c0_X_filteredTest = c0_X_rawTest;


%% Training the ANN %%
Parameters_ANNTraining.type_seeds = 'uniform';%'uniform', 'random' ... anyother idea?
Parameters_ANNTraining.train_sz = 75;
Parameters_ANNTraining.val_sz = 15;
Parameters_ANNTraining.test_sz = 10;
Parameters_ANNTraining.seeds = f_Generate_Seeds(Parameters_ANNTraining.type_seeds, Parameters_ANNTraining.num_seeds);
% Descriptors2Use = []; %all descriptors, otherwise specify, e.g.:{'Grid_5x5_polar_mean', 'Grid_5x5_cartesian_mean', 'Grid_5x5_polar_median', 'Grid_5x5_cartesian_median'};
% Descriptors2Use = {'Grid_5x5_polar_mean', 'Grid_5x5_cartesian_mean', 'Grid_5x5_polar_median', 'Grid_5x5_cartesian_median'};

%% Datasets
%training
ci_X_cip1_training = ci_X_cip1_filteredTrain;
%testing
c0_X_testing = c0_X_filteredTest;
ci_X_cip1_testing = ci_X_cip1_filteredTest;
m_H_ci_testing = m_H_ci_filteredTest;
% cells containing the ANN for Rot and Trans. We used cells in order to
% take advantage of the parallel toolbox
% Centering data
mean_ci_X_cip1_training = mean(ci_X_cip1_training, 2); % rowwise
ci_X_cip1_training = ci_X_cip1_training - mean_ci_X_cip1_training*ones(1, size( ci_X_cip1_training, 2));

ANNs_Translation = cell(num_OF_Datasets_Names, 1);
ANNs_Rotation = cell(num_OF_Datasets_Names, 1);
mean_Descriptors_training = cell(num_OF_Datasets_Names, 1);
for i_method=1:num_OF_Datasets_Names,
    [ANNs_Translation{i_method}, ANNs_Rotation{i_method}, mean_Descriptors_training{i_method}] = f_Train_ANN(ci_X_cip1_training, Descriptors_training.(OF_Datasets_Names{i_method}), Parameters_ANNTraining, Descriptors2Use);
end

%% Predicting with the Trained ANN %%
% here!! extracting/loading descriptors for training
% for now assuming the same!
c0_X_testing_hat = cell(num_OF_Datasets_Names, 1);
ci_X_cip1_testing_hat = cell(num_OF_Datasets_Names, 1);
m_H_testing_hat = cell(num_OF_Datasets_Names, 1);
ANN_Errors = cell(num_OF_Datasets_Names, 1);
% Predicting
for i_method=1:num_OF_Datasets_Names,
    [c0_X_testing_hat{i_method}, ci_X_cip1_testing_hat{i_method}, m_H_testing_hat{i_method}] = f_ANN_Prediction(Descriptors_testing.(OF_Datasets_Names{i_method}), mean_Descriptors_training{i_method},  mean_ci_X_cip1_training, ANNs_Translation{i_method}, ANNs_Rotation{i_method}, Parameters_ANNTraining, Descriptors2Use);
end
%% Plotting
% Computing best/median/wrost trajectories
controlPts = size(c0_X_testing, 2);
for i_method=1:num_OF_Datasets_Names,
    ANN_Errors{i_method} = f_ANN_Plotting_Results(c0_X_testing_hat{i_method}, ci_X_cip1_testing_hat{i_method}, m_H_ci_testing, m_H_testing_hat{i_method}, c0_X_testing, ci_X_cip1_testing, controlPts, Parameters_ANNTraining, Descriptors2Use);
end

%% Final Plots %%
% Computing the Metrics
ind_NN2Plot = 1; % 1 plotting the best seeds, and 2 includes the median, and 3 includes the worst
E_CtrlPts =  cell(ind_NN2Plot, 1);
E_IntMax =  cell(ind_NN2Plot, 1);
E_VTran =  cell(ind_NN2Plot, 1);
E_VRot =  cell(ind_NN2Plot, 1);
numFramesVel = size(ci_X_cip1_testing, 2);
for i_method=1:num_OF_Datasets_Names,
    for i_descriptors=1:length(Descriptors2Use),
        for i_seed = 1:ind_NN2Plot,
            ind2plot = ANN_Errors{i_method}.indices.(Descriptors2Use{i_descriptors})(i_seed);
            [E_CtrlPts{i_seed}, E_IntMax{i_seed}, E_VTran{i_seed}, E_VRot{i_seed}] = f_Compute_Metrics(...
                [], ci_X_cip1_testing_hat{i_method}.(Descriptors2Use{i_descriptors}){ind2plot}, m_H_ci_testing, ci_X_cip1_testing, controlPts);
        end
        Metrics_ControlPts.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_CtrlPts;
        Metrics_IntTrajMax.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_IntMax;
        Metrics_Velocity_Translation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_VTran;
        Metrics_Velocity_Rotation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}) = E_VRot;
    end
end

distances = sqrt(sum(ci_X_cip1_testing(1:3, :).^2));
total_distance = sum(distances);
percentagefactor = (100/total_distance);


% plotting
% 1 - Method/Descriptor Comparison
method_color = 'rbmk';
seed_symbol = 'v+^';
maxInt_symbol = 'hpo';
Descriptors2Use_plotnames = {'CM', 'Cd'};%{'PM', 'CM', 'Pd', 'Cd'};
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
            y = Metrics_IntTrajMax.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}){i_seed};
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
        plot(Metrics_Velocity_Translation.(OF_Datasets_Names{i_method}).(Descriptors2Use{i_descriptors}){bestseed}, 'b-');
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
        h = subplot(1, 2, i_descriptors);
        %h2 = f_Plot_Trajectory(m_H_int, 1:35:700, h, 0.75, 0.5,  'b', 'gt');
        %h2 = f_Plot_Trajectory(m_H_int, 1:35:700, h, 0.75, 0.5,  'b', []);
        h2 = f_Plot_Trajectory_LineOnly(m_H_ci_testing, 1:length(m_H_ci_testing), h, 'b-');
        indmin = ANN_Errors{i_method}.indices.(Descriptors2Use{i_descriptors})(1);
        %f_Plot_Trajectory(m_H_hat{i_method}.(Descriptors2Use{i_descriptors}){indmin}, 1:35:700, h, 0.75, 0.5,  'r', []);
        f_Plot_Trajectory_LineOnly(m_H_testing_hat{i_method}.(Descriptors2Use{i_descriptors}){indmin},...
            1:length(m_H_testing_hat{i_method}.(Descriptors2Use{i_descriptors}){indmin}), h, 'r-');
        xlabel([Descriptors2Use_plotnames{i_descriptors}, ': b=GT, r=est']);
    end
end

%% Boxplots of seeds - ControlPts
figure;
for i_method=1:num_OF_Datasets_Names, %method
    x = zeros(Parameters_ANNTraining.num_seeds, numDescrip);
    for i_descriptors=1:numDescrip,
        x(:, i_descriptors) = percentagefactor*ANN_Errors{i_method}.ControlPts.(Descriptors2Use{i_descriptors});
        ind2rem = x(:, i_descriptors)>100;
        x(ind2rem, i_descriptors) = 100; 
    end
    subplot(2, 2, i_method);
    h = boxplot(x);
    set(gca, 'XTick', 1:numDescrip);
    set(gca, 'XTickLabel', Descriptors2Use_plotnames(1:numDescrip));
    xlabel(OF_Datasets_Names{i_method})
end