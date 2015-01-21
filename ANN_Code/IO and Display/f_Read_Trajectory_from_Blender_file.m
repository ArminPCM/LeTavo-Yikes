%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m_H_c, c0_X, ci_X_cip1, c0_X_int, m_H_int] = f_Read_Trajectory_from_Blender_file(filename, indFrames, Parameters)
%% Function that reads the trajectory data from a text file.
%The file should contain in each row the H ([R|t]) matrix wrt to a global
%reference frame. 
% the output are the centered H, poses with rpy rotations, relative
% motions, integradted rpy poses, and the corresponding H that may differ
% from the original one, but it may be more accurate for comparison wrt the ANN.

%% Nonholonomic trajectory generation
% loading data from the file
nout = max(nargout,1) - 1;
d_camera_motion = load(filename); % in the reference frame of the "data"
numOfFrames = length(indFrames);
d_H_c = cell(numOfFrames, 1);
for i_frame = 1:length(indFrames),
    ind_Frame = indFrames(i_frame);
    d_H_c(i_frame) = {[d_camera_motion(ind_Frame,1:4); ...
        d_camera_motion(ind_Frame,5:8); ...
        d_camera_motion(ind_Frame,9:12); ...
        d_camera_motion(ind_Frame,13:16)]};
end

%% Initial Camera Pose
m_R_c0 = rotox(-pi/2);
m_t_c0 = zeros(3,1);
c0_H_m = [m_R_c0' -m_R_c0'*m_t_c0; 0 0 0 1];
m_H_c0 = inv(c0_H_m);
c0_H_d  = inv(d_H_c{1}); %initial pose, in d reference frame

m_H_c = cell(numOfFrames, 1); % Camera pose in m reference frame
c0_H_c = cell(numOfFrames, 1); % Camera pose in c0 reference frame
%% Center data at origin
c0_X = zeros(6, numOfFrames); % Camera pose in c0 reference frame, rotations given by rpy.
for i_frame = 1:numOfFrames
    % in c0 reference frame
    c0_H_c{i_frame} = c0_H_d*d_H_c{i_frame};
    % In matlab reference frame
    m_H_c{i_frame} = m_H_c0*c0_H_c{i_frame};
    % computing roll, pitch and yaw
    c0_X(1:3, i_frame) = c0_H_c{i_frame} (1:3, 4);
    [r, y, p] = f_R2rpy(c0_H_c{i_frame} (1:3, 1:3));
    c0_X(4:6, i_frame) = mod([r, p, y]+pi, 2*pi)-pi;
end
%% computing Velocities
[ci_X_cip1, m_H_ci] = f_Compute_Velocities(c0_X); %% Plotting - making sure the velocities are correct
% integrating trajectories
[c0_X_int,  m_H_int] = f_Integrated_Trajectory(ci_X_cip1); % integrating trajectories
if Parameters.Visualize,
    h = f_Plot_Trajectory(m_H_c, [1:Parameters.Vis_dx:700], [], 0.75, 1, 'b', 'gt');
%     f_Plot_Trajectory(m_H_int, [1:Parameters.Vis_dx:700], h, 0.75, 1, 'r', 'int');
end