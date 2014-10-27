%% BOX: The Swiss trains example

cities = {'Lausanne',    ...
          'Geneva',      ...
          'Zurich',      ...
          'Neuchatel',   ...
          'Bern'};
      
% Create the train time matrix
D = [ 0   33  128   40   66     
      0    0  158   64  101    
      0    0    0   88   56     
      0    0    0    0   34     
      0    0    0    0    0   ];
  
n  = size(D, 1);
D = D + D';
D = D.^2;

% Do the first variant of classic MDS described in the text (alternatively,
% just call classic_mds)
XX = -1/2 * (D - D(:, 1) * ones(1, n) - ones(n, 1) * D(1, :));
[U, S, V] = svd(XX);
S = S(:, 1:2);
X = sqrt(S')*V';

% Fix one reflection
X(1, :) = -X(1, :);

% Rotate approximately
phi = 130 * pi/180;
R = [ cos(phi) sin(phi)
     -sin(phi) cos(phi) ];
X = R*X;     

% Plot the result
figure(1);
clf;
hold on;
for j = 1:n
    scatter(X(1, j), X(2, j), 'filled');
    text(X(1, j), X(2, j), cities{j});
end
axis equal;
axis tight;
