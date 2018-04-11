f = @(x)power(pi,2)/2*sin(pi*x/2);
p = @(x)1;
q = @(x)power(pi/2,2);
x = [0, 0.1, 0.3, 0.4, 0.777, 0.9, 1]';

x = unique(sort(min(1,max(0, [0; 1; normrnd(0.5, 0.3, 100, 1)]))));
disp(size(x));

[K, f, u] = solve(p, q, f, x);
