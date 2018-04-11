function [K0, f, u] = solve(p, q, f0, x)
    n = size(x,1);
    h = zeros(n);
    for i = uint32(2:n)
        h(i) = x(i)-x(i-1);
    end
    
    K0 = zeros(2*n-1, 2*n-1);
    b = zeros(2*n-1, 1);
    for i = uint32(1:(n-1))
        class(i)
        K = zeros(3, 3);
        f = @(t)dphi_i(i, x, h, t)*dphi_i(i, x, h, t)*p(t)+phi_i(i, x, h, t)*phi_i(i, x, h, t)*q(t);
        if i > 1
            K(1, 1) = integral(f, x(i-1), x(i),'ArrayValued',true);
        end

        integral(@(x)1, x(i), x(i+1), 'ArrayValued', true)
        K(1, 1) = K(1, 1) + integral(f, x(i), x(i+1),'ArrayValued',true);
        
        f = @(t)dphi_i(i, x, h, t)*dphi_i_5(i, x, h, t)*p(t)+phi_i(i, x, h, t)*phi_i_5(i, x, h, t)*q(t);
        K(1, 2) = integral(f, x(i), x(i+1),'ArrayValued',true);
        K(2, 1) = K(1, 2);
        
        f = @(t)dphi_i(i, x, h, t)*dphi_i(i+1, x, h, t)*p(t)+phi_i(i, x, h, t)*phi_i(i+1, x, h, t)*q(t);
        K(1, 3) = integral(f, x(i), x(i+1),'ArrayValued',true);
        K(3, 1) = K(1, 3);
        
        f = @(t)dphi_i_5(i, x, h, t)*dphi_i_5(i, x, h, t)*p(t)+phi_i_5(i, x, h, t)*phi_i_5(i, x, h, t)*q(t);
        K(2, 2) = integral(f, x(i), x(i+1),'ArrayValued',true);
        
        f = @(t)dphi_i(i+1, x, h, t)*dphi_i_5(i, x, h, t)*p(t)+phi_i(i+1, x, h, t)*phi_i_5(i, x, h, t)*q(t);
        K(2, 3) = integral(f, x(i), x(i+1),'ArrayValued',true);
        K(3, 2) = K(2, 3);

        f = @(t)dphi_i(i+1, x, h, t)*dphi_i(i+1, x, h, t)*p(t)+phi_i(i+1, x, h, t)*phi_i(i+1, x, h, t)*q(t);
        K(3, 3) = integral(f, x(i), x(i+1),'ArrayValued',true);
        if i < n-1
            K(3, 3) = K(3, 3) + integral(f, x(i+1), x(i+2),'ArrayValued',true);
        end
        K0((2*i-1):(2*i+1), (2*i-1):(2*i+1)) = K0((2*i-1):(2*i+1), (2*i-1):(2*i+1)) + K;
    end
    
    for i = uint32(1:n)
        f = @(t)f0(t)*phi_i(i, x, h, t);
        if i > 1 
            b(2*i-1) = integral(f, x(i-1), x(i),'ArrayValued',true); 
        end
        if i < n 
            b(2*i-1) = b(2*i-1)+integral(f, x(i), x(i+1),'ArrayValued',true); 
        end
    end
    
    for i = uint32(1:(n-1))
        f = @(t)f0(t)*phi_i_5(i, x, h, t);
        b(2*i) = integral(f, x(i), x(i+1),'ArrayValued',true);
    end
    
    K0
    u = jacobi(K0, b, 1e-5);
end
