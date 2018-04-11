function val = dphi_i(i, x, h, t)
    % Derivative of function phi_i(x)
    n = size(x, 1);   
    if (i > 1) && (x(i-1) <= t) && (t <= x(i))
        val = (3*h(i)+4*t-4*x(i))/power(h(i),2);
    elseif (i < n) && (x(i) < t) && (t <= x(i+1))
        val = (-3*h(i+1)+4*t-4*x(i))/power(h(i+1),2);
    else
        val = 0;
    end
end
