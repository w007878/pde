% 随机初始化一个x
x = max(0.001, min(0.999, normrnd(0.5, 0.25, 100, 1)));
x = unique(sort(x));
n = size(x, 1);

% 问题规格化
f = @(x)pi^2/2*sin(x*pi/2);
p = @(x)1;
q = @(x)pi^2/4;

% 求解
u = linear_element(f, p, q, x);
plot(x, [0; u])
