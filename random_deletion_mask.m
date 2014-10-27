function W = random_deletion_mask(n, n_del)
%%
% W = random_deletion_mask(n, n_del)
%
% Computes the observation mask for an EDM where unobserved distances are
% chosen uniformly at random.
%
% INPUT:  n     ... the number of points in the EDM
%		  n_del ... the number of distances to remove
%
% OUTPUT: W     ... the sought mask
%
% Author: Ivan Dokmanic, 2014

k = randperm(n*(n - 1)/2, n_del);
W = zeros(n);
i = find(tril(ones(n), -1));
W(i(k)) = 1;
W = logical(ones(n) - (W + W'));
