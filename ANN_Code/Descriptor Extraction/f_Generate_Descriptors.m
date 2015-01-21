%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Descriptors = f_Generate_Descriptors(InputParameters, IndSelectedFrames, DescriptorParameters)
%% function that generates the OF descriptors
% Input
% InputParameters: 
% IndSelectedFrames: indices of the optical flow pairs used
% DescriptorParameters: Parameters of the descriptors
%  - DescriptorParameters.EnablingGrid  {0 or 1}: Enables the grid
%          descriptor
%  - DescriptorParameters.GridSize: a 1x2 vector with the grid dimensions
%  - DescriptorParameters.ImageSize: image pixel resolution
%  - DescriptorParameters.GridLimits: 1x4 vector which defines a ROI where
%          the descriptors will be calculated
%  - DescriptorParameters.Units: {'polar' or 'cartesian'} units used for
%           the descriptor. If empty, both cases are considered.
%  - DescriptorParameters.DescriptorFunction: function to calculate the
%           descriptor 'representative' value in each bin based on the
%           patch OF. if empty it estimates both mean and median
%  - DescriptorParameters.DescriptorFunctionName: name (string) for the function,
%           needed for assembling the output structure
%  - DescriptorParameters.MinThreshold: minimum threshold over the 'representative' value
%           of the OF
%  - DescriptorParameters.Normalize: Normalizes each patch based on the
%           patch length
%  - DescriptorParameters.SigmaNoise: (bagging) Adds white Guassian noise
%          with sigma std to each input OF (deprecated, noise is preferably added before passing the OF to this function)
%  - DescriptorParameters.EnablingLumen  {0 or 1}: Enables the lumen-based
%            descriptor.
%  - DescriptorParameters.LumenAngles: Angles for the lumen partitioning (e.g., -pi:pi/2:pi; )
%  - DescriptorParameters.LumenLevels: distance thresholds for the lumen
%            partitioning (e.g., [50 50 50 100 100 Inf];)
% Output
% mxn (specified by the parameters) matrix of m descriptors each of size n

% loading the naem of the OF-variables
if isempty(IndSelectedFrames)
    error('No pairs selected ... pelase select the image pairs to extract descriptors')
end

stringname = ['ProjPairs_', InputParameters.strmethod, '_'];

descriptor_setup_list = f_Setup();

numFrames2Process = length(IndSelectedFrames);
%% Descriptor Calculation
% descriptors_vector = cell(length(descriptor_setup_list), 1);
for i_options=1:length(descriptor_setup_list);
    options = descriptor_setup_list(i_options).options;
    descriptors_i = zeros(2*options.GridSize_1, numFrames2Process);
    dataset = InputParameters.dataset;
    for i_frames=1:numFrames2Process,
        ProjPairs =  f_ReadOF(dataset, stringname, IndSelectedFrames(i_frames));
        desc_i = f_compute_OF_grid_partitioning(ProjPairs.U1, ProjPairs.U2, options, 0);
        descriptors_i(:, i_frames) = desc_i(:);
    end
    Descriptors.(options.name) = descriptors_i;
end

    function options = f_GenerateOptions4GridPartitioning(type)
        % packs the options for the grid partitioning based on the
        % parameters in DescriptorParameters.
        options.type = type;
        options.grid_size = DescriptorParameters.GridSize;
        options.GridSize_1 = prod(DescriptorParameters.GridSize);
        options.image_size = DescriptorParameters.ImageSize;
        options.grid_limits = DescriptorParameters.GridLimits;
        options.units = 'cartesian';
        options.descriptor = @(x)(mean(x, 2));  % default
        options.min_threshold = DescriptorParameters.MinThreshold;
        options.normalize = DescriptorParameters.Normalize;
        options.noise_sigma = DescriptorParameters.SigmaNoise;
        if strcmp(type, 'lumen')
            options.lumen.angles = DescriptorParameters.LumenAngles; %-pi:pi/2:pi;
            options.lumen.levels = DescriptorParameters.LumenLevels; %[50 50 50 100 100 Inf]; % 6 levels
            options.GridSize_1 = length(DescriptorParameters.LumenAngles)*length(DescriptorParameters.LumenLevels)+1;
        end
    end

    function descriptor_setup_list = f_Setup()
        % packs the options for the grid partitioning based on the
        % parameters in DescriptorParameters.
        descriptor_setup_list = [];
        options_vec = [];
        if  DescriptorParameters.EnablingGrid,
            options = f_GenerateOptions4GridPartitioning('grid');
            options.name = ['Grid_', num2str(DescriptorParameters.GridSize(1)), 'x', num2str(DescriptorParameters.GridSize(2)), '_'];
            options_vec(1).options = options;
            %
            if  isempty(DescriptorParameters.Units),
                options_vec(2).options = options;
                options_vec(1).options.units = 'polar';
                options_vec(1).options.name = [options_vec(1).options.name, 'polar_'];
                options_vec(2).options.units = 'cartesian';
                options_vec(2).options.name = [options_vec(2).options.name, 'cartesian_'];
            else
                options_vec(1).options.units = DescriptorParameters.Units;
                options_vec(1).options.name = [options_vec(1).options.name, DescriptorParameters.Units, '_'];
            end
            %
            if  isempty(DescriptorParameters.DescriptorFunction),
                num_opt = length(options_vec);
                options_vec = [options_vec, options_vec];
                for i_opt=1:num_opt,
                    options_vec(i_opt).options.descriptor = @(x)(mean(x, 2));
                    options_vec(i_opt).options.name = [options_vec(i_opt).options.name, 'mean'];
                    options_vec(i_opt+num_opt).options.descriptor = @(x)(median(x, 2));
                    options_vec(i_opt+num_opt).options.name = [options_vec(i_opt+num_opt).options.name, 'median'];
                end
            else
                num_opt = length(options_vec);
                for i_opt=1:num_opt,
                    options_vec(i_opt).options.descriptor = DescriptorParameters.DescriptorFunction;
                    options_vec(i_opt).options.name = [options_vec(i_opt).options.name, DescriptorParameters.DescriptorFunctionName];
                end                
            end
            descriptor_setup_list = options_vec;
        end
        if  DescriptorParameters.EnablingLumen,
            options = f_GenerateOptions4GridPartitioning('lumen');
            options.name = ['Lumen_', num2str(length(DescriptorParameters.LumenAngles)), 'x', num2str(length(DescriptorParameters.LumenLevels)), '+1_'];
            options_vec(1).options = options;
            %
            if  isempty(DescriptorParameters.Units),
                options_vec(2).options = options;
                options_vec(1).options.units = 'polar';
                options_vec(1).options.name = [options_vec(1).options.name, 'polar_'];
                options_vec(2).options.units = 'cartesian';
                options_vec(2).options.name = [options_vec(2).options.name, 'cartesian_'];
            else
                options_vec(1).options.units = DescriptorParameters.Units;
                options_vec(1).options.name = [options_vec(1).options.name, DescriptorParameters.Units, '_'];
            end
            %
            if  isempty(DescriptorParameters.DescriptorFunction),
                num_opt = length(options_vec);
                options_vec = [options_vec, options_vec];
                for i_opt=1:num_opt,
                    options_vec(i_opt).options.descriptor = @(x)(mean(x, 2));
                    options_vec(i_opt).options.name = [options_vec(i_opt).options.name, 'mean'];
                    options_vec(i_opt+num_opt).options.descriptor = @(x)(median(x, 2));
                    options_vec(i_opt+num_opt).options.name = [options_vec(i_opt+num_opt).options.name, 'median'];
                end
            else
                num_opt = length(options_vec);
                for i_opt=1:num_opt,
                    options_vec(i_opt).options.descriptor = DescriptorParameters.DescriptorFunction;
                    options_vec(i_opt).options.name = [options_vec(i_opt).options.name, DescriptorParameters.DescriptorFunctionName];
                end
            end
            descriptor_setup_list = [descriptor_setup_list options_vec];
        end
    end    
% end of main function 
end

function ProjPairs =  f_ReadOF(dataset, stringname, index)
        name = [stringname, num2str(index), '_', num2str(index+1)];
        load(dataset, name);
        eval(['ProjPairs = ', name, ';']);
    end