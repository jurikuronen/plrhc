function [Ghathc, t] = hc(data, Ghator, g, penalty, direction)
% Apply the second phase global hill-climbing on the OR solution Ghator computed using hcor.m.
% Juri Kuronen (2018)
%
% Input:
% - data: discrete Nxd matrix
% - Ghator: logical dxd matrix, the OR solution Ghathc computed using hcor.m
% - g: gamma parameter for extended BIC (0 - classical BIC)
% - penalty: 2x1 information vector; first element is penalty type
%   (0 = no Lq-penalty, 1 = L1-penalty, 2 = L2-penalty), second element is lambda value
% - direction: 0 - start from OR graph (recommended), 1 - start from empty graph
%
% Output:
% - Ghathc: logical dxd matrix, the second phase HC solution computed from Ghator
% - t: total time taken

    t = tic;
    var_noc = max(data);
    d = size(data, 2);
    if direction == 1; Ghathc = sparse(d, d); else; Ghathc = Ghator; end

    % Get current Markov blankets and compute their scores
    mbs = arrayfun(@(i) find(Ghathc(:, i))', 1:d, 'uniformoutput', false);
    mb_scores = arrayfun(@(i) calc_ebicscore(data, i, mbs{i}, [], g, penalty, var_noc), 1:d); 

    % Compute 'gain' of adding or removing each edge. Store this information in matrix form as
    % [i j mbscore_new(i) mbscore_new(j) gain]
    edges = [];
    for i = 1:d-1
        for j = find(Ghathc(i, i+1:end)) + i
            [mbi, mbj] = get_new_mbs(i, j, mbs{i}, mbs{j});
            mbscorei = calc_ebicscore(data, i, mbi, [], g, penalty, var_noc);
            mbscorej = calc_ebicscore(data, j, mbj, [], g, penalty, var_noc);
            edges = [edges; [i j mbscorei mbscorej mbscorei + mbscorej - mb_scores(i) - mb_scores(j)]];
        end
    end

    if isempty(edges); Ghathc = sparse(d, d); t = toc(t); return; end
    [max_gain, max_idx] = max(edges(:, 5));
    while max_gain > 0
        % Add or remove best edge
        i = edges(max_idx, 1); j = edges(max_idx, 2);
        mbscorei = mb_scores(i); mbscorej = mb_scores(j);                   % Save old scores
        mb_scores(i) = edges(max_idx, 3); mb_scores(j) = edges(max_idx, 4); % Update scores
        edges(max_idx, 3) = mbscorei; edges(max_idx, 4) = mbscorej;         % Resetting this change would give back old scores
        edges(max_idx, 5) = -edges(max_idx, 5);                             % New gain is naturally the inverse
        [mbi, mbj] = get_new_mbs(i, j, mbs{i}, mbs{j});                     % Update mbs after the change
        mbs{i} = mbi; mbs{j} = mbj;

        % Apply corrections to gains for neighbors of i or j
        edges = apply_changes(data, edges, i, j, 1, mbs, mb_scores, g, penalty, var_noc);
        edges = apply_changes(data, edges, i, j, 2, mbs, mb_scores, g, penalty, var_noc);
        edges = apply_changes(data, edges, j, i, 1, mbs, mb_scores, g, penalty, var_noc);
        edges = apply_changes(data, edges, j, i, 2, mbs, mb_scores, g, penalty, var_noc);

        [max_gain, max_idx] = max(edges(:, 5));
    end

    Ghathc = sparse(d, d);
    for i = 1:d
        Ghathc(mbs{i}, i) = 1;
        Ghathc(i, mbs{i}) = 1;
    end
    t = toc(t);

end

function edges = apply_changes(data, edges, node, node2, pos, mbs, mb_scores, g, penalty, var_noc)
% After node-node2 edge was modified by main algorithm, apply corrections to neighbors of node and node2.

    for ind = find(edges(:, pos) == node)'
        neighbor = edges(ind, mod(pos, 2) + 1);
        if neighbor == node2; continue; end
        [mb, ~] = get_new_mbs(node, neighbor, mbs{node}, mbs{neighbor});
        edges(ind, pos + 2) = calc_ebicscore(data, node, mb, [], g, penalty, var_noc);
        edges(ind, 5) = sum(edges(ind, 3:4)) - mb_scores(node) - mb_scores(neighbor);
    end

end

function [mbi mbj] = get_new_mbs(i, j, mbi, mbj)

    if sum(mbi == j) > 0
        mbi = setdiff(mbi, j);
        if nargout == 2; mbj = setdiff(mbj, i); end
    else
        mbi = [mbi j];
        if nargout == 2; mbj = [mbj i]; end
    end

end
