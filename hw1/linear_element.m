function [u] = linear_element(f, p, q, x_s)
    x = @(i)x_s(i + 1);
    
    n = size(x_s, 1) - 1;
    h = zeros(n, 1);
    for i = 1:n
        h(i) = x(i) - x(i - 1);
    end

    xi = @(i, t)(t - x(i - 1))/h(i);
    
    K = zeros(n, n);
    u = zeros(n + 1, 1);
    for i = 2:n
        func = @(t)(p(x(i-1) + h(i)*xi(i, t))/h(i) + h(i)*q(x(i-1)+h(i)*xi(i, t))*power((1-xi(i, t)),2));
        K(i - 1, i - 1) = K(i - 1, i - 1) + integral(func, 0, 1);

        func = @(t)(p(x(i-1) + h(i)*xi(i, t))/h(i) + h(i)*q(x(i-1)+h(i)*xi(i, t))*power(xi(i, t),2));
        K(i, i) = K(i, i) + integral(func, 0, 1);

        func = @(t)(p(x(i-1) + h(i)*xi(i, t))/h(i) + h(i)*q(x(i-1)+h(i)*xi(i, t))*(xi(i, t) .* (1-xi(i, t))));
        
        K(i, i - 1) = K(i, i - 1) + integral(func, 0, 1);
        K(i - 1, i) = K(i - 1, i) + integral(func, 0, 1);
    end
    
    b = zeros(n, 1);
    for i = 1:n
        func = @(t)(f(x(i-1) + h(i)*xi(i, t)) .* xi(i, t));
        b(i) = h(i) * integral(func, 0, 1);
    end
    for i = 2:n
        func = @(t)(f(x(i-1) + h(i)*xi(i, t)) .* (1 - xi(i, t)));
        b(i - 1) = b(i - 1) + h(i) * integral(func, 0, 1);
    end
    
    u = pinv(K) * b;
end