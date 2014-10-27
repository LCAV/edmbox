MAX_ITER = 50;
n = 20;
d = 2;

tic;
for i = 1:MAX_ITER	
	X = rand(d, n);
	D = edm(X);
	[~, D_est] = alternating_descent(D, d);
end
t = toc
