%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figHandle = f_Plot_Trajectory(m_H_c, indFrames, figHandle, param_cam_scale, param_axis_scale,  param_axis_color, param_axis_preffix)

if isempty(figHandle)
    figHandle = figure;
end
% figure(figHandle);
if ~ishold, 
    hold on;
end
for i=indFrames,
    f_3Dcamera(m_H_c{i} , param_axis_color ,param_cam_scale, 2); % pyramidal case
    if ~isempty(param_axis_preffix),
        f_3Dframe(m_H_c{i}, param_axis_color, param_axis_scale, [param_axis_preffix, num2str(i)]); %axis
    else
        plot3(m_H_c{i}(1, 4), m_H_c{i}(2, 4), m_H_c{i}(3, 4), ['o' param_axis_color]);
    end
end
axis equal;