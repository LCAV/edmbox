function X = classic_mds(D, dim)
%%
% X = classic_mds(D, dim)
%
% Reconstructs the (centered) point set from a noisy distance matrix.
%
% INPUT:   D   ... measured Euclidean Distance Matrix (EDM) (n by n)
%          dim ... target embedding dimension
%
% OUTPUT:  X   ... (dim by n) list of point coordinates
%
% Author: Ivan Dokmanic, 2014

n = size(D, 1);
I = eye(n);
J = I - 1/n*ones(n);

[U, S, V] = svd(-J*D*J/2);
S = S(1:dim, :);
X = sqrt(S)*V';
