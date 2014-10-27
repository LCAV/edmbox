function D = sdr_complete_D(t_D, W, dim)
%%
% D = sdr_complete_D(t_D, W, dim)
%
% Completes a partially observed D (with observed entries assumed noiseless)
% using the semidefinite relaxation (SDR).
%
% INPUT:  t_D ... input partial EDM
%         W   ... observation mask (1 for observed, 0 for non-observed entries)
%         dim ... desired embedding dimensions
%
% OUTPUT: D   ... completed EDM
%
% Author: Ivan Dokmanic, 2014

n = size(t_D, 1);

% The old SVD code to do the following is also good.
x = -1/(n + sqrt(n));
y = -1/sqrt(n);
V = [y*ones(1, n-1); x*ones(n-1) + eye(n-1)];
e = ones(n, 1);

cvx_quiet('true');
cvx_begin sdp
    variable G(n-1, n-1) symmetric;
    maximize trace(G);
    B = V*G*V';
    E = diag(B)*e' + e*diag(B)' - 2*B;
    subject to
        E .* W == t_D .* W;
        G >= 0;
cvx_end
cvx_quiet('false');

% Can I just use what's above inside the cvx block? Don't do that, do rank
% thresholding here!
B = V*G*V';
D = diag(B)*e' + e*diag(B)' - 2*B;

