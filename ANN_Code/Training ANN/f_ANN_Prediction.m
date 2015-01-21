%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c0_X_hat, ci_X_cip1_hat, m_H_hat] = f_ANN_Prediction(Testing_Descriptors, Descriptors_mean, ci_X_cip1_mean, ANN_Translation, ANN_Rotation, Parameters, Descriptors2Use, paralel_option)

if nargin < 7 || isempty(Descriptors2Use),
    Descriptors2Use = fieldnames(Testing_Descriptors);
end

display('Predicting Motion!')
%% Prediccion
c0_X_hat_i  =  cell(Parameters.num_seeds, 1);
m_H_hat_i  =  cell(Parameters.num_seeds, 1);
for i_descriptor = 1:length(Descriptors2Use),
    display(Descriptors2Use{i_descriptor});
    X_testing = Testing_Descriptors.(Descriptors2Use{i_descriptor}) - Descriptors_mean.(Descriptors2Use{i_descriptor})*ones(1, size(Testing_Descriptors.(Descriptors2Use{i_descriptor}), 2));
    ANN_Translation_i = ANN_Translation.(Descriptors2Use{i_descriptor});
    ANN_Rotation_i      = ANN_Rotation.(Descriptors2Use{i_descriptor});
    if paralel_option,
        parfor i_seed = 1:Parameters.num_seeds,
            Rel_Translation_i = sim(ANN_Translation_i{i_seed}, X_testing);
            Rel_Rotation_i      = sim(ANN_Rotation_i{i_seed},      X_testing);
            ci_X_cip1_hat_i(i_seed) = {ci_X_cip1_mean*ones(1, size(Rel_Translation_i, 2)) + [Rel_Translation_i; Rel_Rotation_i]};
            [c0_X_int,  m_H_int_i] = f_Integrated_Trajectory(ci_X_cip1_hat_i{i_seed}); % integrating trajectories
            c0_X_hat_i(i_seed) = {c0_X_int};
            m_H_hat_i(i_seed) = {m_H_int_i};
        end
    else
        for i_seed = 1:Parameters.num_seeds,
            Rel_Translation_i = sim(ANN_Translation_i{i_seed}, X_testing);
            Rel_Rotation_i      = sim(ANN_Rotation_i{i_seed},      X_testing);
            ci_X_cip1_hat_i(i_seed) = {ci_X_cip1_mean*ones(1, size(Rel_Translation_i, 2)) + [Rel_Translation_i; Rel_Rotation_i]};
            [c0_X_int,  m_H_int_i] = f_Integrated_Trajectory(ci_X_cip1_hat_i{i_seed}); % integrating trajectories
            c0_X_hat_i(i_seed) = {c0_X_int};
            m_H_hat_i(i_seed) = {m_H_int_i};
        end
    end
    ci_X_cip1_hat.(Descriptors2Use{i_descriptor}) = ci_X_cip1_hat_i;
    c0_X_hat.(Descriptors2Use{i_descriptor}) = c0_X_hat_i;
    m_H_hat.(Descriptors2Use{i_descriptor}) = m_H_hat_i;
end