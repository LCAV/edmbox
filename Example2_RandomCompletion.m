%% Compares the performance of several algorithms for EDM completion, when
%  the unobserved entries are randomly chosen in the matrix.

n = 20; % Number of points in the set
d = 2;  % Embedding dimension

% Minimal and maximal number of deletions
n_del_min  = 20;
n_del_step = 5;
n_del_max  = 150;
n_del      = n_del_min:n_del_step:n_del_max;

% Number of random point sets per deletion level
n_config = 500;

% Number of random masks for every point configuration
n_mask = 1;

% Run the simulation
methods = {'Alternating Descent', ...
           'Rank Alternation', ...
           'Semidefinite Relaxation'};
err = zeros(3, numel(n_del));
success = zeros(3, numel(n_del));
SUCC_TOL = 1e-2;

parfor i_del = 1:numel(n_del)
    
    err_in = zeros(3, 1);
    success_in = zeros(3, 1);
    for i_config = 1:n_config
        
        fprintf('#(Deletions) %d in the range %d-%d, configuration %d/%d\n', n_del(i_del), n_del_min, n_del_max, i_config, n_config);

        X = rand(d, n);      % Point set
        G = X'*X;            % Gramian
        D = edm(X, X); % EDM
        
        for i_mask = 1:n_mask
            W = random_deletion_mask(n, n_del(i_del));
            
            % Alternating Coordinate Descent
            [~, E] = alternating_descent(D .* W, d);
            err_in(1) = err_in(1) + norm(E - D, 'fro');
            success_in(1) = success_in(1) + (norm(E - D, 'fro') < SUCC_TOL*norm(D, 'fro'));
            
            % Rank EDM Complete
            E = rank_complete_edm(D, W, d, 0);
            
            err_in(2) = err_in(2) + norm(E - D, 'fro');
            success_in(2) = success_in(2) + (norm(E - D, 'fro') < SUCC_TOL*norm(D, 'fro'));
            
            % Semidefinite Relaxation 1
            E = sdr_complete_edm(D, W, d);
            err_in(3) = err_in(3) + norm(E - D, 'fro');
            success_in(3) = success_in(3) + (norm(E - D, 'fro') < SUCC_TOL*norm(D, 'fro'));
        end
    end
    err(:, i_del) = err_in;
    success(:, i_del) = success_in;
end

err = err / n_config / n_mask;
success = success / n_config / n_mask;

%% Plotting

figure(1);
clf;

plot(n_del, success', 'LineWidth', 2);
ylabel('Success percentage');
xlabel('Number of deletions');
legend(methods, 'Location','SouthWest');

axis tight;
grid on;

