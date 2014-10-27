%% Compares the performance of several algorithms on Metric Multidimensional 
%  Unfolding, for various numbers of acoustic events.

m = 20; % Number of microphones
d = 3;  % Embedding dimension

% Minimal and maximal number of acoustic events
n_events_min =  5;
n_events_step = 1;
n_events_max = 30;
n_events = n_events_min:n_events_step:n_events_max;

% Number of random point sets per event number
n_config = 100;

% Run the simulation
methods = {'A Method by Crocco', ...
           'Alternating Descent', ...
           'Rank alternation', ...
           'Semidefinite Relaxation'};
err = zeros(4, numel(n_events));
success = zeros(4, numel(n_events));
SUCC_TOL = 1e-2;

parfor i_events = 1:numel(n_events)
    
    k = n_events(i_events);
    
    err_in = zeros(4, 1);
    success_in = zeros(4, 1);
    for i_config = 1:n_config
        
        fprintf('#(Events) %d in the range %d-%d, configuration %d/%d\n', k, n_events_min, n_events_max, i_config, n_config);

        X = rand(d, m);        % Microphones
        Y = rand(d, k);        % Acoustic events
        D = edm([X Y], [X Y]); % EDM
        W = mdu_mask(m, k);

        % Crocco et al.
        [~, ~, E] = unfold_crocco(sqrt(edm(X, Y)), 1);
        err_in(1) = err_in(1) + norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro');
        success_in(1) = success_in(1) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < SUCC_TOL*norm(D(1:m, 1:m), 'fro'));

        % Alternating Coordinate Descent
        [~, E] = alternating_descent(D .* W, d);
        err_in(2) = err_in(2) + norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro');
        success_in(2) = success_in(2) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < SUCC_TOL*norm(D(1:m, 1:m), 'fro'));

        % Rank EDM Complete
        E = rank_complete_edm(D, W, d, 0);

        err_in(3) = err_in(3) + norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro');
        success_in(3) = success_in(3) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < SUCC_TOL*norm(D(1:m, 1:m), 'fro'));

        % Semidefinite Relaxation 1
        E = sdr_complete_edm(D, W, d);
        err_in(4) = err_in(4) + norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro');
        success_in(4) = success_in(4) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < SUCC_TOL*norm(D(1:m, 1:m), 'fro'));
    end
    err(:, i_events) = err_in;
    success(:, i_events) = success_in;
end

err = err / n_config;
success = success / n_config;

%% Plotting

figure(1);

plot(n_events, success', 'LineWidth', 2);
ylabel('Success percentage');
xlabel('Number of acoustic events');
legend(methods, 'Location','East');

axis tight;
grid on;

