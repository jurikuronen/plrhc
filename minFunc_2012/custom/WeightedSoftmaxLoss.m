function [nll, g, H] = WeightedSoftmaxLoss(w, X, y, k, weights)

[n, p] = size(X);
w = reshape(w, [p k-1]);
w(:, k) = zeros(p, 1);

Z = sum(exp(X * w), 2);

nll = -sum(weights .* (sum(X .* w(:, y).', 2) - log(Z)));

if nargout > 1
    g = zeros(p, k-1);
    for c = 1:k-1
        g(:, c) = -sum(bsxfun(@times, X .* repmat((y == c) - exp(X * w(:, c)) ./ Z, [1 p]), weights));
    end
    g = reshape(g, [p * (k - 1) 1]);
end

if nargout > 2
    H = zeros(p * (k - 1));
    SM = exp(X * w(:, 1:k-1)) ./ repmat(Z, [1 (k - 1)]);
    for c1 = 1:k-1
        for c2 = 1:k-1
            D = weights .* SM(:, c1) .* ((c1 == c2) - SM(:, c2));
            H(((p * (c1 - 1) + 1)):(p * c1), ((p * (c2 - 1) + 1)):(p * c2)) = X' * diag(sparse(D)) * X;
        end
    end
end
