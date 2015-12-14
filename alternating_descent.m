function [X, D] = alternating_descent(t_D, dim)
%%
% [X, D] = alternating_descent(t_D, dim)
%
% INPUT:  t_D ... (n by n) noisy and incomplete observation of an EDM for n points
%         dim ... (scalar) target embedding dimension
%
% OUTPUT: X   ... (dim by n) the reconstructed point set
%         D   ... (n by n)   the nearest EDM with embedding dimension dim
%         
% Authors: Reza Parhizkar and Ivan Dokmanic, 2014


MAX_ITER = 50;
[n, m] = size(t_D);
if (n ~= m)
    error('The inpur matrix t_D needs to be a square matrix!')
end

L = eye(n) - 1/n*ones(n,1)*ones(n,1).';

X0 = zeros(n, dim);
X = X0;

% extract connectivity information

connected_nodes = cell(1, n);
connectivity_vector = cell(1, n);
sensor_ind_list = 1:n;

for ind = 1 : n
    sensor_ind = sensor_ind_list(ind);
    connected_nodes{ind} = find(t_D(sensor_ind, :) ~= 0);
    connectivity_vector{ind} = t_D(sensor_ind, connected_nodes{ind});
    connectivity_vector{ind} = connectivity_vector{ind}(:);
end

for iter_ind = 1 : MAX_ITER
    for ind = 1 : n
        
        X_connected = X(connected_nodes{ind}, :);
        for coor_ind = 1 : dim
            a = 4 * size(X_connected, 1);
            b = 12 * sum((X(ind,coor_ind) -  X_connected(:,coor_ind)));
            c = 4 * sum(sum((repmat(X(ind,:), size(X_connected, 1), 1) - X_connected).^2, 2) + 2 * (X(ind,coor_ind) -  X_connected(:,coor_ind)).^2 - connectivity_vector{ind});
            d = 4 * sum( (X(ind,coor_ind) -  X_connected(:,coor_ind)) .* (sum((repmat(X(ind,:), size(X_connected, 1), 1)- X_connected).^2, 2) - connectivity_vector{ind}));
            
            DeltaX_vec = cubicfcnroots(a,b,c,d);
            DeltaX_vec = real(DeltaX_vec(abs(imag(DeltaX_vec)) < 1e-15));
            
            if isempty(DeltaX_vec)
                continue;
            end
            
            cc = zeros(1, length(DeltaX_vec));
            for ii = 1 : length(DeltaX_vec)
                cc(ii) = sum(((X(ind,coor_ind)+DeltaX_vec(ii)-X_connected(:,coor_ind)).^2 ...
                       + sum((repmat(X(ind,:), size(X_connected, 1), 1) - X_connected).^2, 2) ... 
                       - (X(ind, coor_ind) - X_connected(:, coor_ind)).^2 ...
                       - connectivity_vector{ind}).^2);
            end
            [~, cc_min_ind] = min(cc);
            DeltaX = DeltaX_vec(cc_min_ind);
            X(ind,coor_ind) = X(ind,coor_ind) + DeltaX;
            
            X = L * X;
        end
    end
end

X = L * X;
[k, kk] = find(triu(ones(n),1));
D = zeros(n);
D(k + n * (kk - 1)) = sum((X(k,:) - X(kk,:)).^2, 2);
D(kk + n * (k - 1)) = D(k + n * (kk - 1));

