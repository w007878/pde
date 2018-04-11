function val = phi_i_5(i, x, h, t)
    % function phi_i.5(x) by interpolation
    n = size(x, 1);
    if i < n && x(i) <= t && t <= x(i+1)
        val = 4*(t-x(i))/h(i+1) * (1-(t-x(i))/h(i+1));
    else
        val = 0;
    end
end
