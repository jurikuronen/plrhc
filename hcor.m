function [Ghator t] = hcor(data, Gstar, g, penalty, depth, parallel)
% Apply the first phase local hill-climbing algorithm for each variable to
% learn the OR solution.
% Juri Kuronen (2018)
%
% Input:
% - data: discrete Nxd matrix
% - Gstar: logical dxd matrix, the graph G^star with edges E^star computed
%   using plr.m (optional, can be left as [])
% - g: gamma parameter for extended BIC (0 - classical BIC)
% - penalty: 2x1 information vector; first element is penalty type
% - depth: the search space will contain all nodes within graph-theoretic-distance 'depth'
%   from a given node with respect to E^star (only used if Gstar is provided)
% - parallel: with positive value, uses MATLAB's parallelization
%
% Output:
% - Ghator: logical dxd matrix, the OR solution
% - t: time taken per node or total time taken if computed in parallel

    d = size(data, 2);
    var_noc = max(data);
    plr = (size(Gstar, 1) == d);
    Ghator = sparse(d, d);

    if parallel > 0
        t = tic;
        parfor i = 1:d
            if plr; searchspace = compute_searchspace(Gstar, i, depth); else; searchspace = []; end;
            [mb, ~] = learn_mb(data, i, searchspace, g, penalty, var_noc);
            tmp = zeros(d, 1); tmp(mb) = 1;
            Ghator(:, i) = tmp;
        end
        t = toc(t);
    else
        t = zeros(d, 1);
        for i = 1:d
            ti = tic;
            if plr; searchspace = compute_searchspace(Gstar, i, depth); else; searchspace = []; end;
            [mb, ~] = learn_mb(data, i, searchspace, g, penalty, var_noc);
            Ghator(mb, i) = 1;
            t(i) = toc(ti);
        end
    end

    Ghator = ((Ghator ~= 0) + (Ghator' ~= 0)) > 0;

end
