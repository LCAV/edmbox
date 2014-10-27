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

L = eye(n) - 1/n*ones(n);

% Add the possibility of a warm start (e.g. via classic MDS)
X0 = zeros(n, dim);
D = zeros(n);
X = X0;


for iter_ind = 1 : MAX_ITER
    connectivity_vector = [];
    connected_nodes = [];
    sensor_ind_list = 1:n;

    for ind = 1 : n
        
        % find sensors that are connected to sensor i (non-zero elements of
        % t_D)
        sensor_ind = sensor_ind_list(ind);
        connected_nodes = find(t_D(sensor_ind, :) ~= 0);
        connectivity_vector = t_D(sensor_ind, connected_nodes);
        connectivity_vector = connectivity_vector(:);
        X_connected = X(connected_nodes, :);
        
        for coord = 1 : dim
            
            a = 4 * size(X_connected, 1);
            b = 12 * sum((X(ind,coord) -  X_connected(:,coord)));
            c = 4 * sum(sum((repmat(X(ind,:), size(X_connected, 1), 1)- X_connected).^2, 2) + 2 * (X(ind,coord) -  X_connected(:,coord)).^2 - connectivity_vector);
            d = 4 * sum( (X(ind,coord) -  X_connected(:,coord)) .* (sum((repmat(X(ind,:), size(X_connected, 1), 1)- X_connected).^2, 2) - connectivity_vector));
            
            DeltaX_vec = cubicfcnroots(a,b,c,d);
            DeltaX_vec = real(DeltaX_vec(abs(real(DeltaX_vec)- DeltaX_vec) < 1e-15));
            DeltaX_vec = [DeltaX_vec];
            
            cc = [];
            for ii = 1 : length(DeltaX_vec)
                cc(ii) = sum(((X(ind,coord)+DeltaX_vec(ii)-X_connected(:,coord)).^2 + sum((repmat(X(ind,:), size(X_connected, 1), 1)- X_connected).^2, 2) - (X(ind, coord) - X_connected(:, coord)).^2 - connectivity_vector).^2);
            end
            [~, cc_min_ind] = min(cc);
            DeltaX = DeltaX_vec(cc_min_ind);
            X(ind,coord) = X(ind,coord) + DeltaX;
            X = L * X;
        end
    end
end
X = L * X;
[k,kk] = find(triu(ones(n),1));
D = zeros(n);
D(k + n * (kk - 1)) = sum((X(k,:) - X(kk,:)).^2, 2);
D(kk + n * (k - 1)) = D(k + n * (kk - 1));d