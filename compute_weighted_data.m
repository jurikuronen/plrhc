function [data, weights] = compute_weighted_data(data, node, mb, noc_max)
% compute_weighted_data.m
% Computes a weighted data set from common outcome permutations.
%
% Input:
% - data: discrete Nxd matrix
% - node: index of the node of interest
% - mb: indices of the Markov blanket nodes
% - noc_max: maximum outcome value
%
% Output: 
% - data: data set with weighted rows
% - weights: weight of each row

    % In case of problems with MEX, avoid using a weighted data set
    % by swapping the code commentings below.
    [data, weights] = mex_compute_weighted_data(data, node, mb, noc_max);
    %weights = ones(size(data, 1), 1);

end
