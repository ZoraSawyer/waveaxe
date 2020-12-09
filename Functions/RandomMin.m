function [ minX ] = RandomMin( input_set )

% MATLAB function RandomMin:

%   Input variables: 
%           input_set: a vector containig the input set of numbers

%   Output variables: 
%           minX: the minimum value of the randomly selected subset

subset_size = 6;                                    % size of the randomly selected subset (6 for this example)
index = randperm(length(input_set),subset_size);    % generate random indices
subset = input_set(index);                          % select the random subset
minX = min(subset);                                 % min values of the subset

end