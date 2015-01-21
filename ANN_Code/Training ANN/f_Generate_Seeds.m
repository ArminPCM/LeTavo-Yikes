%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G. Puerto-Souza
%  gustavhafen@gmail.com
%  Astra Lab
%
%  Updated: Jan 20th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seeds = f_Generate_Seeds(type_seeds, num_seeds)
% function that generates the seeds for the neural networks.  
% Possible improvements: more types. Bounds for interval, etc

if strcmp(type_seeds, 'uniform'),
    seeds = 1/(num_seeds+1):1/(num_seeds+1):num_seeds/(num_seeds+1);
end

if strcmp(type_seeds, 'random'),
    seeds = rand(num_seeds,1);
end

end