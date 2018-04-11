function val = phi_i(i, x, h, t)
    n = size(x, 1);
    if i > 1 && x(i-1) <= t && t <= x(i)
        val = (2*(x(i)-t)/h(i)-1) * ((-t+x(i))/h(i)-1);
    elseif i < n && x(i) < t && t <= x(i+1)
        val = (2*(-x(i)+t)/h(i+1)-1) * ((t-x(i))/h(i+1)-1);        
    else
        val = 0;
    end
end
