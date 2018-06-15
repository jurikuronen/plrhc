function [Gstar t] = plr(data, g, penalty, parallel)
% Use penalized likelihood ratio tests to form edges E^star,
% which contains the candidate neighborhoods.
% Juri Kuronen (2018)
%
% Input:
% - data: discrete Nxd matrix
% - g: gamma parameter for extended BIC (0 - classical BIC)
% - penalty: 2x1 information vector; first element is penalty type
%   (0 = no Lq-penalty, 1 = L1-penalty, 2 = L2-penalty), second element is lambda value
% - parallel: with positive value, uses MATLAB's parallelization
%
% Output:
% - Gstar: logical dxd matrix, the graph G^star with edges E^star
% - t: total time taken

    t = tic;
    var_noc = max(data);
    d = size(data, 2);
    Gstar = sparse(d, d);

    empty_mb_scores = arrayfun(@(i) calc_ebicscore(data, i, [], [], g, penalty, var_noc), 1:d); 

    if parallel > 0
        parfor i = 1:d
            tmp = zeros(d, 1);
            for j = i+1:d
                if calc_ebicscore(data, i, j, [], g, penalty, var_noc) > empty_mb_scores(i) || calc_ebicscore(data, j, i, [], g, penalty, var_noc) > empty_mb_scores(j)
                    tmp(j) = 1;
                end
            end
            Gstar(i, :) = tmp;
        end
    else
        for i = 1:d
            for j = i+1:d
                if calc_ebicscore(data, i, j, [], g, penalty, var_noc) > empty_mb_scores(i) || calc_ebicscore(data, j, i, [], g, penalty, var_noc) > empty_mb_scores(j)
                    Gstar(i, j) = 1;
                end
            end
        end
    end

    Gstar = ((Gstar ~= 0) + (Gstar' ~= 0)) > 0;
    t = toc(t);

end
