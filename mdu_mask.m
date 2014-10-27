function W = mdu_mask(m, k)
%%
% W = mdu_mask(m, k)
%
% Creates the masking matrix for multidimensional unfolding to be used with
% semidefinite programs. We use the calibration terminology (microphones and
% acoustic events).
%
% INPUT:  m ... number of microphones
%         k ... number of acoustic events (calibration events)
%
% OUTPUT: W ... the sought mask
%
% Author: Ivan Dokmanic, 2014

W = [zeros(m)   ones(m, k) ;
     ones(k, m) zeros(k)  ];
W = logical(W);