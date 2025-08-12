function [U, S, V, H_proj] = projectToLowRankHankel(H, r)
% Low-rank approximation and Hankel projection of matrix H
% Inputs:
%   H : Input Hankel matrix
%   r : Desired rank for low-rank approximation
% Outputs:
%   U, S, V : SVD components of original H (truncated)
%   H_proj : Hankel-projected matrix with rank-r approximation

    % === Step 1: SVD and rank truncation ===
    [U, S, V] = svd(H, 'econ');
    S(r+1:end, r+1:end) = 0;
    H_approx = U * S * V';

    % === Step 2: Projection onto Hankel structure (anti-diagonal averaging) ===
    [m, n] = size(H_approx);
    H_proj = zeros(m, n);

    for k = 2:(m + n)
        [i, j] = find(bsxfun(@plus, (1:m)', (1:n)) == k);
        avg = mean(H_approx(sub2ind([m, n], i, j)));
        H_proj(sub2ind([m, n], i, j)) = avg;
    end
end
