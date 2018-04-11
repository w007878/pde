function x = jacobi(A, b, eps)
    D = diag(diag(A));
    U = triu(A) - D;
    L = tril(A) - D;
    x = ones(size(A, 1), 1);
    x0 = x;
    x = pinv(D)*(L+U)*x+pinv(D)*b;
    while norm(x - x0) > eps
        x0 = x;
        x = pinv(D)*(L+U)*x+pinv(D)*b;
    end
end
