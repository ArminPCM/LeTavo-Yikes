%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figHandle = f_Plot_Trajectory_LineOnly(m_H_c, indFrames, figHandle,colorSpec)

if isempty(figHandle)
    figHandle = figure;
end
% figure(figHandle);
if ~ishold, 
    hold on;
end
for i=1:length(indFrames),
    trans(:,i) = m_H_c{indFrames(i)}(1:3,4);
end

plot3(trans(1,:),trans(2,:),trans(3,:),colorSpec, 'LineWidth', 2)


axis equal;