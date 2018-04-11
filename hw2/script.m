% Initialize the problem, use the same one to homework1
f = @(x)power(pi,2)/2*sin(pi*x/2);
p = @(x)1;
q = @(x)power(pi/2,2);

% Random choose x by normal distribution.
x = unique(sort(min(1,max(0, [0; 1; normrnd(0.5, 0.3, 100, 1)]))));
n = size(x,1);
h = zeros(n);
    
% compute the difference between two point
    
for i = uint32(2:n)
    h(i) = x(i)-x(i-1);
end

% Solve the problem
[K, f, u] = solve(p, q, f, x, h);

y = plotans(x, h, u);
plot(x(2:(n-1)), y);
