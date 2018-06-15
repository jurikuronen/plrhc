function [mbhat, mbhat_score] = learn_mb(data, node, searchspace, g, penalty, var_noc)
% Apply the local hill-climbing algorithm to learn the Markov blanket of a given node.
% Juri Kuronen (2018)
%
% Input:
% - data: discrete Nxd matrix
% - node: column index of the node
% - searchspace: constrained search space for the given node (optional, can be left as [])
% - g: gamma parameter for extended BIC (0 - classical BIC)
% - penalty: 2x1 information vector; first element is penalty type
%   (0 = no Lq-penalty, 1 = L1-penalty, 2 = L2-penalty), second element is lambda value
% - var_noc: maximum outcome values of all variables
%
% Output:
% - mbhat: column indices of the learned Markov blanket
% - mbhat_score: eBIC score for the learned Markov blanket

    mbhat = [];
    if var_noc(node) == 1; mbhat_score = -inf; return; end
    d = size(data, 2);
    if length(searchspace) == 0; searchspace = setdiff(1:d, node)'; end

    % Keep track of earlier beta estimates to speed-up L-BFGS algorithm
    B = cell(d + 1, var_noc(node) - 1);

    % Score of an empty Markov blanket
    [mbhat_score beta_hat] = calc_ebicscore(data, node, [], [], g, penalty, var_noc);

    for c = 1:length(beta_hat)
        B{1, c} = beta_hat(c); % Intercept
    end

    mbhat_loc = [];
    for i = 1:length(searchspace)
        x = searchspace(i);
        [score, beta_hat] = calc_ebicscore(data, node, x, [], g, penalty, var_noc);
        for c = 1:size(beta_hat, 2)
            B{1, c} = beta_hat(1, c);
            B{x + 1, c} = beta_hat(2:end, c);
        end
        if score > mbhat_score
            mbhat_score = score;
            mbhat = x;
            mbhat_loc = i;
        end
    end
    if length(mbhat) == 0; return; end;
    searchspace(mbhat_loc) = [];

    CONT = 1;
    while CONT == 1
        CONT = 0;
        [mbcand_topscore, mbcand_toploc, B] = find_best_addition(data, node, searchspace, mbhat, B, g, penalty, var_noc);
        if mbcand_topscore > mbhat_score
            mbhat = [mbhat searchspace(mbcand_toploc)];
            searchspace(mbcand_toploc) = [];
            mbhat_score = mbcand_topscore;
            CONT = 1;
        end

        DEL = 1;
        while DEL == 1 && CONT == 1 && length(mbhat) > 2
            DEL = 0;
            [mbcand_topscore, mbcand_toploc, B] = find_best_deletion(data, node, mbhat, B, g, penalty, var_noc);
            if mbcand_topscore > mbhat_score
                mbhat_score = mbcand_topscore;
                searchspace(end + 1) = mbhat(mbcand_toploc);
                mbhat(mbcand_toploc) = [];
                DEL = 1;
            end
        end
    end
end

function [mbcand_topscore, mbcand_toploc, B] = find_best_addition(data, node, searchspace, mb, B, g, penalty, var_noc)
    mbcand_topscore = -inf;
    mbcand_toploc = -1;
    for i = 1:length(searchspace)
        mbcand = [mb searchspace(i)];
        [mbcand_score beta_hat] = calc_ebicscore(data, node, mbcand, cell2mat(B([1 mbcand + 1], :)), g, penalty, var_noc);
        % Update estimates.
        for c = 1:size(beta_hat, 2)
            B{1, c} = beta_hat(1, c); % Intercept
            v = 2;
            for x = mbcand
                v2 = var_noc(x) - 2;
                B{x + 1, c} = beta_hat(v:v+v2, c);
                v = v + v2 + 1;
            end
        end
        if mbcand_score > mbcand_topscore
            mbcand_topscore = mbcand_score;
            mbcand_toploc = i;
        end
    end
end

function [mbcand_topscore, mbcand_toploc, B] = find_best_deletion(data, node, mb, B, g, penalty, var_noc) 
    mbcand_topscore = -inf;
    mbcand_toploc = -1;
    for i = 1:length(mb)
        mbcand = setdiff(mb, mb(i));
        [mbcand_score beta_hat] = calc_ebicscore(data, node, mbcand, cell2mat(B([1 mbcand + 1], :)), g, penalty, var_noc);
        % Update estimates.
        for c = 1:size(beta_hat, 2)
            B{1, c} = beta_hat(1, c); % Intercept.
            v = 2;
            for x = mbcand
                v2 = var_noc(x) - 2;
                B{x + 1, c} = beta_hat(v:v+v2, c);
                v = v + v2 + 1;
            end
        end
        if mbcand_score > mbcand_topscore
            mbcand_topscore = mbcand_score;
            mbcand_toploc = i;
        end
    end
end

