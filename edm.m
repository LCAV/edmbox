function D = edm(X, Y)
%%
% D = edm(X, Y)
% 
% Computes the EDM between the point lists X and Y.
%
% INPUT:  X, Y ... (dim by m) and (dim by n) point lists
%
% OUTPUT: D    ... (m by n) Euclidean distance matrix for X and Y
%
% Note: D = edm(X) is the same as D = edm(X, X)
%
% Author: Ivan Dokmanic, 2014


if nargin < 2
	Y = X;
end

norm_X2 = sum(X.^2); 
norm_Y2 = sum(Y.^2); 

D = bsxfun(@plus, norm_X2', norm_Y2) - 2*X'*Y;
