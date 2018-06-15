function [score beta_hat] = calc_ebicscore(data, node, mb, beta_init, g, penalty, var_noc)
% Calculate the score of a node given a Markov blanket with a logistic regression model.
% Juri Kuronen (2018)
%
% Input:
% - data: discrete Nxd matrix
% - node: column index of the node
% - mb: column indices of the Markov blanket nodes
% - beta_init: initial guess for beta-parameters (optional, can be left as [])
% - g: gamma parameter for extended BIC (0 - classical BIC)
% - penalty: 2x1 information vector; first element is penalty type
%   (0 = no Lq-penalty, 1 = L1-penalty, 2 = L2-penalty), second element is lambda value
% - var_noc: maximum outcome values of all variables
%
% Output:
% - score: eBIC score of the node given its Markov blanket
% - beta_hat: beta parameter maximum likelihood estimates

    node_noc = var_noc(node);
    if node_noc <= 1; score = -inf; beta_hat = beta_init; return; end
    [N, d] = size(data);

    % Use a weighted data set to speed up computations.
    [weighted_data, weights] = compute_weighted_data(data, node, mb, max(var_noc([node mb])));
    y = weighted_data(:, node);
    X = ones(size(weights, 1), 1);
    if length(mb) > 0
        for x = 1:length(mb)
            for l = 2:var_noc(mb(x))
                X = [X weighted_data(:, mb(x)) == l];
            end
        end
    end

    % Handle wrong beta_init input.
    if size(beta_init, 1) ~= size(X, 2); beta_init = zeros(size(X, 2), node_noc - 1); end

    % Loss function for logistic regression
    mfOptions.Display = 'off';
    mfOptions.Method = 'lbfgs';
    if node_noc > 2 % Multinomial LR
        beta_init = beta_init(:);
        funObj = @(beta)WeightedSoftmaxLoss(beta, X, y, node_noc, weights);
    else
        y = sign((y == y(1)) - 0.5); % Map y to {-1, 1}.
        funObj = @(beta)WeightedLogisticLoss(beta, X, y, weights);
    end

    % Lq-penalty
    if penalty(1) == 0 % Don't use Lq-penalty
        [beta_hat ML] = minFunc(funObj, beta_init, mfOptions);
    elseif penalty(1) == 1 % L1-penalty
        mfOptions.verbose = 0;
        lambda = penalty(2) * [zeros(1, node_noc - 1); ones(size(X, 2) - 1, node_noc - 1)];
        lambda = lambda(:);
        [beta_hat ML] = L1General2_PSSgb(funObj, beta_init, lambda, mfOptions);
    else % L2-penalty
        lambda = penalty(2) * [zeros(1, node_noc - 1); ones(size(X, 2) - 1, node_noc - 1)];
        lambda = lambda(:);
        funObj2 = @(beta)penalizedL2(beta, funObj, lambda);
        [beta_hat ML] = minFunc(funObj2, beta_init, mfOptions);
    end

    if node_noc > 2; beta_hat = reshape(beta_hat, size(X, 2), node_noc - 1); end

    % eBIC score
    score = -sum(ML) - (node_noc - 1) * size(X, 2) * log(N) / 2 - length(mb) * g * log(d - 1);

end
