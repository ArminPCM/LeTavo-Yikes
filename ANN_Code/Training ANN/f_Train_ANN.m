%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ANN_Translation, ANN_Rotation, Descriptors_mean] = f_Train_ANN(ci_X_cip1, Descriptors, Parameters, Descriptors2Use, paralel_option)
%% Function that trains the ANN for the relative motions (ci_X_cip1) and descriptors specified in Descriptors2Use (otherwise for all the fields available in Descriptors) and based on the Parameters.
% Inputs:
% ci_X_cip1: 6xn matrix of relative motions
% Descriptors: mxn matrix of descriptors
% Parameters: parameters for the ANN
% Descriptors2Use: fieldnames of the descriptors used to estimate the ANN
%
%Outputs:
% ANN struct containig n fields (one per descriptor in Descriptors2Use),
% each field contains a cell with Parameters.num_seeds ANN obtained
% corresponding to the initial weigths specified in Parameters.seeds.

if nargin<4 || isempty(Descriptors2Use),
    Descriptors2Use = fieldnames(Descriptors);
end

num_descriptors = length(Descriptors2Use);
% Training data Y
Y_Translation = ci_X_cip1(1:3, :);
Y_Rotation = ci_X_cip1(4:6, :);
% ANNs
ANN_Translation_i = cell(Parameters.num_seeds, 1);
ANN_Rotation_i = cell(Parameters.num_seeds, 1);
seeds = Parameters.seeds;
train_sz = Parameters.train_sz;
val_sz = Parameters.val_sz;
test_sz = Parameters.test_sz;
for i_descriptor=1:num_descriptors,
    display(Descriptors2Use{i_descriptor});
    % removing NAN columns!
    num_rows = size(Descriptors.(Descriptors2Use{i_descriptor}), 1);
    Descriptors_mean_i = NaN*zeros(num_rows, 1);
    for i_mean=1:num_rows,
        ind_nan = ~isnan(Descriptors.(Descriptors2Use{i_descriptor})(i_mean, :));
        Descriptors_mean_i(i_mean) = mean(Descriptors.(Descriptors2Use{i_descriptor})(i_mean, ind_nan));
    end
    X_training = Descriptors.(Descriptors2Use{i_descriptor}) - Descriptors_mean_i*ones(1, size(Descriptors.(Descriptors2Use{i_descriptor}), 2)); % training data X
    if paralel_option
        parfor i_seed = 1:Parameters.num_seeds,
            display(i_seed);
            ANN_Translation_i(i_seed) = {train_NNs(X_training, Y_Translation, train_sz, val_sz, test_sz, seeds(i_seed))};
            ANN_Rotation_i(i_seed)      = {train_NNs(X_training, Y_Rotation, train_sz, val_sz, test_sz, seeds(i_seed))};
        end
    else
        for i_seed = 1:Parameters.num_seeds,
            display(i_seed);
            ANN_Translation_i(i_seed) = {train_NNs(X_training, Y_Translation, train_sz, val_sz, test_sz, seeds(i_seed))};
            ANN_Rotation_i(i_seed)      = {train_NNs(X_training, Y_Rotation, train_sz, val_sz, test_sz, seeds(i_seed))};
        end
    end
    ANN_Translation.(Descriptors2Use{i_descriptor}) = ANN_Translation_i;
    ANN_Rotation.(Descriptors2Use{i_descriptor}) = ANN_Rotation_i;
    Descriptors_mean.(Descriptors2Use{i_descriptor}) = Descriptors_mean_i;
end