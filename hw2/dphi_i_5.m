function val = dphi_i_5(i, x, h, t)
    % Derivative of function phi_i.5(x)
    n = size(x, 1);
    if i < n && x(i) <= t && t <= x(i+1)
        val = 4*(h(i+1)-2*t+2*x(i))/power(h(i+1),2);
    else
        val = 0;
    end
end
