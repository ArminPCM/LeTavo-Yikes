%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OF_partitioning = f_compute_OF_grid_partitioning(XGrid, YGrid, options, plot_option)
% Types of OF paritioning (options)
% options.type = 'grid' or 'lumen'
% options.image_size = [xl, yl] - grid of axb
% options.grid_size = [c, d] - grid of axb
% options.grid_limits = [a, b] - grid applied on interval [c,d]
% options.units = 'polar' or 'cartesian'
% options.descriptor = pointer to function
% options.min_threshold - minimum OF value to be considered (values smaller than these are set to 0)


%% 0) setting parameters
if nargin < 4,
    plot_option = 0;
end

[xl, yl] = deal(options.image_size(1), options.image_size(2));
[image_size_x, image_size_y] = deal(options.image_size(1), options.image_size(2));
[xmin, xmax, ymin, ymax] = deal(options.grid_limits(1), options.grid_limits(2), options.grid_limits(3), options.grid_limits(4));
%% 1) computing Optical Flow
num_pts = size(YGrid, 2);
YGrid = YGrid + options.noise_sigma*randn(2, num_pts);
XGrid = XGrid + options.noise_sigma*randn(2, num_pts);
%
ind2remY = (YGrid(1,:) < 0)|(YGrid(1,:)>xl)|(YGrid(2,:) < 0)|(YGrid(2,:)>yl);
ind2remX = (XGrid(1,:) < 0)|(XGrid(1,:)>xl)|(XGrid(2,:) < 0)|(XGrid(2,:)>yl);
ind2rem = logical(ind2remX+ind2remY);
YGrid(:, ind2rem) = [];
XGrid(:, ind2rem) = [];
num_pts = size(YGrid, 2);
%
OF_vectors = YGrid - XGrid;
r_vector = zeros(num_pts,1);
theta_vector = zeros(num_pts,1);
for i =1:num_pts,
    r_vector(i) = norm(OF_vectors(:, i));
    theta_vector(i) = atan2(OF_vectors(2, i), OF_vectors(1, i));
end
ind_set_2_0 = r_vector < options.min_threshold;
OF_vectors(:, ind_set_2_0) = zeros(2, sum(ind_set_2_0));
r_vector(ind_set_2_0) = 0;
theta_vector(ind_set_2_0) = 0;

% Grid type
if strcmp(options.type, 'grid'),
    [grid_size_x, grid_size_y] = deal(options.grid_size(1), options.grid_size(2));
% %     OF_partitioning = zeros(2, grid_size_x*grid_size_y);
    OF_partitioning = NaN*zeros(2, grid_size_x*grid_size_y);
    % computing the grid positions
    dx =(xmax-xmin)/grid_size_x;
    dy = (ymax-ymin)/grid_size_y;
    x_limits = unique([xmin:dx:xmax]);
    y_limits = unique([ymin:dy:ymax]);
    % last index includes the leftover at the border of the image
    x_limits(end:end+1) = [xmax image_size_x] ;
    y_limits(end:end+1) = [ymax image_size_y] ;
    % computing respective bins for the OF
    ind_x = zeros(num_pts,1);
    ind_y = zeros(num_pts,1);
    for i = 1:num_pts,
        ind_x(i) = find(XGrid(1,i) <= x_limits,1)-1;
        ind_y(i) = find(XGrid(2,i) <= y_limits,1)-1;
    end        
    %computing bins
    if plot_option,
        figure,
        [aa, bb] = meshgrid(x_limits, y_limits);
        mesh(aa, bb, zeros(7,7)); hold on; view(0, 90); axis ij;
        %         plot(XGrid(1, :), XGrid(2, :), 'b+');
        %         plot(YGrid(1, :), YGrid(2, :), 'r+')
        line([XGrid(1, 1:end); XGrid(1, 1:end)+OF_vectors(1, 1:end)], [XGrid(2, 1:end); XGrid(2, 1:end)+OF_vectors(2, 1:end)], 'Color', 'r');
    end
    
    for i=1:grid_size_x,
        for j=1:grid_size_y,
            ind_ij = ((ind_x == i)&(ind_y == j));
            % resampling OF based on the descriptor function
            if sum(ind_ij)>0,
                %%               plot(XGrid(1, ind_ij), XGrid(2, ind_ij), 'rs')
                if strcmp(options.units, 'cartesian')
                    OF_partitioning(:, (i-1)*grid_size_y+j) = options.descriptor(OF_vectors(:, ind_ij)); % make sure the function is applied column-wise!
                end
                if strcmp(options.units, 'polar')
                    r_mean = mean(r_vector(ind_ij));
                    % checking the correct angle:
                    theta_mean1 = mean(mod(theta_vector(ind_ij), 2*pi)); % from 0 - 2pi
                    theta_mean2 = mean(mod(theta_vector(ind_ij)-pi/2, 2*pi)+pi/2); % from pi/2 - pi/2
                    theta_mean3 = mean(mod(theta_vector(ind_ij)+pi, 2*pi)-pi);  % from -pi to pi
                    theta_mean4 = mean(mod(theta_vector(ind_ij)+pi/2, 2*pi)-pi/2); % from -pi/2 - pi/2
                    theta_mean = mode([theta_mean1, theta_mean2, theta_mean3, theta_mean4]);
                    theta_mean = mod(theta_mean+pi/2, 2*pi)-pi/2;
                    OF_partitioning(: , (i-1)*grid_size_y+j) = [r_mean; theta_mean];
                end
            end
        end
    end
    if plot_option,
        [XX, YY] = meshgrid(xmin+dx/2:dx:xmax, ymin+dy/2:dy:ymax );
        XY  = [XX(:) YY(:)];
        OF_partitioning2 = OF_partitioning(:);
        if strcmp(options.units, 'cartesian')
            line([XY(:, 1) XY(:, 1)+OF_partitioning2(1:2:end)]', [XY(:, 2) XY(:, 2)+OF_partitioning2(2:2:end)]', 'Color', 'b');
        end
        if strcmp(options.units, 'polar')
            line([XY(:, 1) XY(:, 1)+OF_partitioning2(1:2:end).*cos(OF_partitioning2(2:2:end))]', [XY(:, 2) XY(:, 2)+OF_partitioning2(1:2:end).*sin(OF_partitioning2(2:2:end))]', 'Color', 'b');
        end
    end
    %% Normalization for grid 
    if options.normalize,
        if strcmp(options.units, 'polar'),
            % computing maximum length in square: 
            norm_r = sqrt(dx*dx+dy*dy);
            norm_t = 2*pi;
            OF_partitioning(2, :) = mod(OF_partitioning(2, :)+pi/2, 2*pi)-pi/2; %from 0 to 2pi
            OF_partitioning = [1/norm_r 0; 0 1/norm_t]*OF_partitioning;
        end
        if strcmp(options.units, 'cartesian'),
            % normalizing x
            OF_partitioning = [1/dx 0; 0 1/dy]*OF_partitioning;

        end
    end
end

if strcmp(options.type, 'lumen'),
    %% parameters for lumen-based partitioning
    % options.lumen.levels = offset for the different length levels (different radius for the grid)
    % options.lumen.angles = number of directions (equally spaced angles from -pi to pi)
    % options.lumen.center = data with the lumen centers
    % options.lumen.radius = data with the lumen radius
    num_lumen_levels = length(options.lumen.levels);
    num_lumen_angles = length(options.lumen.angles);
% %     OF_partitioning = zeros(2, 1+num_lumen_levels*(num_lumen_angles-1));
	OF_partitioning = NaN*zeros(2, 1+num_lumen_levels*(num_lumen_angles-1));%%%
    % computing the lumen grid
    lumen_center = options.lumen.center;
    lumen_radius = options.lumen.radius;
    
    dlevels = cumsum([0 lumen_radius options.lumen.levels]);
    dangles = options.lumen.angles;
    
    % computing respective bins for the OF
    lumen_difference = XGrid - lumen_center*ones(1, num_pts);
    
    r_vectors_lumen = zeros(num_pts, 1);
    theta_vector_lumen = zeros(num_pts, 1);
    ind_r_lumen = zeros(num_pts, 1);
    ind_t_lumen  = zeros(num_pts, 1);
    for i = 1:num_pts,
        r_vectors_lumen(i) = norm(lumen_difference(:, i));
        theta_vector_lumen(i) =atan2(lumen_difference(2, i), lumen_difference(1, i));
        ind_r_lumen(i) =  find(r_vectors_lumen(i) <= dlevels,1);
        ind_t_lumen(i) = find(theta_vector_lumen(i)<= dangles,1);
    end          
    %
    % adjusting indices
    ind_r_lumen = ind_r_lumen - 1;
    ind_t_lumen = ind_t_lumen - 1;
    %computing bins
    if plot_option,
        figure,
        angles = -pi:pi/18:pi;
        circle_1 = [cos(angles); sin(angles)];
        plot(lumen_radius*circle_1(1, :)+lumen_center(1), lumen_radius*circle_1(2, :)+lumen_center(2), 'k-'); %First class
        axis equal; axis ij; axis([0 xl, 0, yl]); hold on;
        for i_l=2:num_lumen_levels
            plot(dlevels(i_l+1)*circle_1(1, :)+lumen_center(1), dlevels(i_l+1)*circle_1(2, :)+lumen_center(2), 'b-'); %First class
        end
        Xline = [0 xl; 0 0];
        for i_l=1:num_lumen_angles-1
            Rot = [cos(dangles(i_l)) -sin(dangles(i_l)); sin(dangles(i_l)) cos(dangles(i_l))];
            XlineR = Rot*Xline+lumen_center*ones(1, size(Xline,2));
            plot(XlineR(1,:), XlineR(2, :), 'r-');
        end
        plot(XGrid(1, :), XGrid(2, :), 'b+');
        plot(YGrid(1, :), YGrid(2, :), 'r+')
        line([XGrid(1, 1:end); XGrid(1, 1:end)+OF_vectors(1, 1:end)], [XGrid(2, 1:end); XGrid(2, 1:end)+OF_vectors(2, 1:end)], 'Color', 'r');
    end
    
    %% Case 1, the lumen!
    ind_ij = ind_r_lumen == 1;
    % resampling OF based on the descriptor function
    if sum(ind_ij)>0, % optical flow in the lumen sector!
        %%               plot(XGrid(1, ind_ij), XGrid(2, ind_ij), 'rs')
        if strcmp(options.units, 'cartesian')
            OF_partitioning(:, 1) = options.descriptor(OF_vectors(:, ind_ij)); % make sure the function is applied column-wise!
        end
        if strcmp(options.units, 'polar')
            r_mean = mean(r_vector(ind_ij));
            % checking the correct angle:
            theta_mean1 = mean(mod(theta_vector(ind_ij), 2*pi)); % from 0 - 2pi
            theta_mean2 = mean(mod(theta_vector(ind_ij)-pi/2, 2*pi)+pi/2); % from pi/2 - pi/2
            theta_mean3 = mean(mod(theta_vector(ind_ij)+pi, 2*pi)-pi);  % from -pi to pi
            theta_mean4 = mean(mod(theta_vector(ind_ij)+pi/2, 2*pi)-pi/2); % from -pi/2 - pi/2
            theta_mean = mode([theta_mean1, theta_mean2, theta_mean3, theta_mean4]);
            OF_partitioning(: , 1) = [r_mean; theta_mean];
        end
    end
    %% Other cases
    for i=1:num_lumen_levels, %
        for j=1:num_lumen_angles-1, %note we are assuming that angles are given from [-pi, pi], including pi! change later if needed...
            ind_ij = ((ind_r_lumen == i+1)&(ind_t_lumen == j));
            % resampling OF based on the descriptor function
            if sum(ind_ij)>0,
                %%               plot(XGrid(1, ind_ij), XGrid(2, ind_ij), 'rs')
                if strcmp(options.units, 'cartesian')
                    OF_partitioning(:, 1+(i-1)*(num_lumen_angles-1)+j) = options.descriptor(OF_vectors(:, ind_ij)); % make sure the function is applied column-wise!
                end
                if strcmp(options.units, 'polar')
                    r_mean = mean(r_vector(ind_ij));
                    % checking the correct angle:
                    theta_mean1 = mean(mod(theta_vector(ind_ij), 2*pi)); % from 0 - 2pi
                    theta_mean2 = mean(mod(theta_vector(ind_ij)-pi/2, 2*pi)+pi/2); % from pi/2 - pi/2
                    theta_mean3 = mean(mod(theta_vector(ind_ij)+pi, 2*pi)-pi);  % from -pi to pi
                    theta_mean4 = mean(mod(theta_vector(ind_ij)+pi/2, 2*pi)-pi/2); % from -pi/2 - pi/2
                    theta_mean = mode([theta_mean1, theta_mean2, theta_mean3, theta_mean4]);
                    OF_partitioning(: , 1+(i-1)*(num_lumen_angles-1)+j) = [r_mean; theta_mean];
                end
            end
        end
    end
        
    if plot_option,
        XY  = lumen_center';        
        dlevels2 = dlevels(2:end); % removing 0 and the first radius
        dlevels2(isinf(dlevels2)) = dlevels(end-1)+100; % assigning the maximum image length in case of infinite margins. 
        middle_pts = (dlevels2(2:end) + dlevels2(1:end-1))/2;
        middle_angles = (dangles(2:end) + dangles(1:end-1))/2;
        Xline = [middle_pts; 0*middle_pts];        
        
        for i_l=1:length(middle_pts),
            point_l = Xline(:, i_l);
            for j_l=1:length(middle_angles)
                Rot = [cos(middle_angles(j_l)) -sin(middle_angles(j_l)); sin(middle_angles(j_l)) cos(middle_angles(j_l))];
                XY = [XY; (Rot*point_l+lumen_center)'];                
            end                        
        end
        plot(XY(:, 1), XY(:,2), 'kd');
        OF_partitioning2 = OF_partitioning(:);
        if strcmp(options.units, 'cartesian')
            line([XY(:, 1) XY(:, 1)+OF_partitioning2(1:2:end)]', [XY(:, 2) XY(:, 2)+OF_partitioning2(2:2:end)]', 'Color', 'b');
        end
        if strcmp(options.units, 'polar')
            line([XY(:, 1) XY(:, 1)+OF_partitioning2(1:2:end).*cos(OF_partitioning2(2:2:end))]', [XY(:, 2) XY(:, 2)+OF_partitioning2(1:2:end).*sin(OF_partitioning2(2:2:end))]', 'Color', 'b');
        end
    end
    
        %% Normalization
    if options.normalize, 
        if strcmp(options.units, 'polar'),
            % computing the maximum distance from the center to image
            % borders
            max_dist = max([xmax - lumen_center(1), lumen_center(1), ymax - lumen_center(2), lumen_center(2)]);
            dlevels2 = dlevels(2:end); % removing 0 and the first radius
            dlevels2(dlevels2 > max_dist) = [];
            dlevels2(end+1) = max_dist;
            norm_r = max(dlevels2(2:end) - dlevels2(1:end-1));
            norm_t = 2*pi;
            OF_partitioning(2, :) = mod(OF_partitioning(2, :)+pi/2, 2*pi)-pi/2; %from 0 to 2pi
            OF_partitioning = [1/norm_r 0; 0 1/norm_t]*OF_partitioning;
        end
        if strcmp(options.units, 'cartesian'),
            % computing the maximum distance from the center to image
            % borders
            max_dist = max([xmax - lumen_center(1), lumen_center(1), ymax - lumen_center(2), lumen_center(2)]);
            dlevels2 = dlevels(2:end); % removing 0 and the first radius
            dlevels2(dlevels2 > max_dist) = [];
            dlevels2(end+1) = max_dist;
            norm_r = max(dlevels2(2:end) - dlevels2(1:end-1));
            OF_partitioning = [1/norm_r 0; 0 1/norm_r]*OF_partitioning;
        end

    end
end


