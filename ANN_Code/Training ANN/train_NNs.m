function net = train_NNs(inputs, targets, train_sz, val_sz, test_sz, rand_seed_val)

% Solve an Input-Output Fitting problem with a Neural Network
% Script generated by NFTOOL
% Created Wed Sep 04 14:08:50 CDT 2013
%
% This script creates neural networks based on the inputs and targets
% given;
%
% Inputs: data which serves as inputs to the neural network
% Targets: data which serves as outputs from the neural network
% Train_sz: percent of the training set which will be used as training data
% Val_sz: percent of the training data which will be used as validation
% Test_sz: percent of the training data which will be used as test

% Process inputs
if(nargin==2)
    train_sz = 80;
    val_sz = 15;
    test_sz = 5;
    rand_seed_val = 0; %Means that you will have random weights (i.e., default behavior)
end

if(nargin==5)
   rand_seed_val=0;
end

% Get apppropriate layer sizes
[num_input_lay_nodes, ~] = size(inputs);
[num_output_lay_nodes,~] = size(targets);


% Create a Fitting Network
num_hidden_lay_nodes = 3; %usually 2n+1
net = fitnet(num_hidden_lay_nodes);
% net = newrb(inputs, targets, 0.0, 0.01, num_hidden_lay_nodes);
% Choose Input and Output Pre/Post-Processing Functions
% For a list of all processing functions type: help nnprocess
% net.inputs{1}.processFcns = {'mapminmax'};
% net.outputs{2}.processFcns = {'mapminmax'};
net.inputs{1}.processFcns =  {'fixunknowns', 'removeconstantrows', 'processpca'}; 
net.outputs{2}.processFcns = {'fixunknowns', 'removeconstantrows', 'processpca'};
% net.inputs{1}.processFcns =  {'fixunknowns', 'removeconstantrows'}; 
% net.outputs{2}.processFcns = {'fixunknowns', 'removeconstantrows'};
% net.inputs{1}.processFcns =  {'fixunknowns'}; 
% net.outputs{2}.processFcns = {'fixunknowns'};
%Set configuration values to initialize the layers of the neural network to
%ready for training
net = configure(net, inputs, targets);

%Set the seed to something we know.
% rand_seed_val
% rng(rand_seed_val);
% rng(rand_seed_val*(2^32-1));
rng(rand_seed_val*(2^10-1));

%Set "random" weights
num_input_lay_nodes = net.inputs{1}.processedSize; %% TAVO!!
% hidden_weights=0.1*randn(num_hidden_lay_nodes, num_input_lay_nodes);
% hidden_b_weights=0.1*randn(num_hidden_lay_nodes,1);

hidden_weights=0.1*randn(num_hidden_lay_nodes, num_input_lay_nodes);
hidden_b_weights=0.1*randn(num_hidden_lay_nodes,1);

net.IW{1,1} = hidden_weights;
net.b{1,1} = hidden_b_weights;

% output_weights=0.1*randn(num_output_lay_nodes, num_hidden_lay_nodes);
% output_b_weights=0.1*randn(num_output_lay_nodes,1);
% net.LW{2,1} = output_weights;
% net.b{2,1} = output_b_weights;


% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'divideblock';%'divideint';  % Divide data in an interleaved fashion
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = train_sz/100;
net.divideParam.valRatio = val_sz/100;
net.divideParam.testRatio = test_sz/100;

% For help on training function 'trainlm' type: help trainlm
% For a list of all training functions type: help nntrain
net.trainFcn = 'trainlm';  % Levenberg-Marquardt for paper comparison
%net.trainFcn = 'trainscg'; %with a quickness or low memory
%net.trainFcn = 'trainbfg';

%Set the transfer functions of the neural network
%Note: turning all transfer functions into purelin makes two matrices
%(i.e., traditional least squares regression only)
% Choose {logsig,purelin,tansig} <- for others "help nntransfer" but these
% are most common.
% Default is {logsig, purelin}.
net.layers{1}.transferFcn='logsig';%'purelin';
net.layers{2}.transferFcn='tansig';

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
% Choose{'mae', 'sae', 'sse', 'mse'}; where {'m':mean; 's':sum; 'a':absolute}
net.performFcn = 'mse';  % Mean squared error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
  'plotregression', 'plotfit'};

%Set training parameters of the network
net.trainParam.max_fail=50; %number of validation failures for early stopping
net.trainParam.goal=1e-24;     %lowest value of the error you set above
%net.trainParam.mu = 0.001; %Establishes an initial relation between GD and GN methods


% Train the Network
net.trainParam.showWindow=0;
net.trainParam.epochs = 5000;
[net,tr] = train(net,inputs,targets);

% Test the Network
% outputs = net(inputs);

% Don't want to take time to compute this just now.
% errors = gsubtract(targets,outputs);
% performance = perform(net,targets,outputs);

% Recalculate Training, Validation and Test Performance
% trainTargets = targets .* tr.trainMask{1};
% valTargets = targets  .* tr.valMask{1};
% testTargets = targets  .* tr.testMask{1};

% Don't want to take the time to process this
% trainPerformance = perform(net,trainTargets,outputs);
% valPerformance = perform(net,valTargets,outputs);
% testPerformance = perform(net,testTargets,outputs);

% View the Network
% view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, plotfit(net,inputs,targets)
%figure, plotregression(targets,outputs)
%figure, ploterrhist(errors)

end
