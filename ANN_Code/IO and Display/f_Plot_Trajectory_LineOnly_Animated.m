%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figHandle = f_Plot_Trajectory_LineOnly_Animated(m_H_c, indFrames, figHandle,colorSpec)

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

X = trans(1,1);
Y = trans(2,1);
Z = trans(3,1);
h = plot3(X, Y, Z,colorSpec);
axis equal;
for i=2:length(trans(1,:)),
    X = trans(1,1:i);
    Y = trans(2,1:i);
    Z = trans(3,1:i);
    set(h, 'XData', X, 'YData', Y, 'ZData', Z);
%     drawnow;
    pause(0.01);
end

