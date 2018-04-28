f = @(p)10;
fd = @(p)sqrt(sum(p.^2,2))-1;
huniform = @(x, y)ones(size(x, 1), 1);
[p, t] = distmesh2d(fd, huniform, 0.1,[-1,-1;1,1],[]);
[uh, in] = poisson(f, fd, 0.5, p, t);

[p, t] = recmesh2d(0.1, [-1,-1;1,1]);
[uh, in] = poisson(f, fd, 0.5, p, t);
