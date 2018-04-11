function y = plotans(x, h, u)
    % get y-coordidate of the plot
    n = size(x, 1);
    y = zeros(n - 2, 1);
    for i = 2:(n-1)
        y(i-1) = u(2*i-3)*phi_i(i-1,x,h,x(i))+u(2*i-2)*phi_i_5(i-1,x,h,x(i));
        y(i-1) = y(i-1)+u(2*i-1)*phi_i(i,x,h,x(i))+u(2*i)*phi_i_5(i,x,h,x(i));
        y(i-1) = y(i-1)+u(2*i+1)*phi_i(i+1,x,h,x(i));
    end
end