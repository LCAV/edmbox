function D = rank_complete_edm(t_D, W, dim, gram)

%%
% D = sdr_complete_D(t_D, W, dim)
%
% Completes a partially observed D (with observed entries assumed noiseless)
% by alternating between enforcing the rank and enforcing the observed
% entries.
%
% INPUT:  t_D  ... input partial EDM
%         W    ... observation mask (1 for observed, 0 for non-observed entries)
%         dim  ... desired embedding dimensions
%         gram ... rank-threshold the Gram matrix, not the EDM
%
% OUTPUT: D    ... completed EDM
%
% Author: Ivan Dokmanic, 2014


% Stopping criteria
MAX_ITER = 1000;
MAX_TOL  = 1e-7;


W      = logical(W);
n      = size(t_D, 1);
mean_d = mean(mean(sqrt(t_D)));
D      = t_D;
D(W)   = mean_d;
I      = logical(eye(n));
J      = I - (1/n)*ones(n);

while true
    for i = 1:MAX_ITER
        % Save D to monitor convergence
        D_old = D;

        % Do rank thresholding
        if gram == 1
            G = -1/2 * J * D * J;
            [U, S, V] = svd(G);
            S(dim+1:end, dim+1:end) = 0;
            G = U * S * V';
            D = diag(G)*ones(1, n) + ones(n, 1)*diag(G)' - 2*G;
        else
            [U, S, V] = svd(D);
            S(dim+3:end, dim+3:end) = 0;
            D = U * S * V';
        end

        change1 = norm(D_old - D, 'fro');

        % Enforce observed entries, zero diagonal and positivity
        D(W) = t_D(W);
        D(I) = 0;
        D(D<0) = 0;

        change2 = norm(D_old - D, 'fro');

        if (change1 < MAX_TOL) && (change2 < MAX_TOL)
            break;
        end
    end
    
    if gram < 2
        break;
    else
        gram = 1;
    end
end
