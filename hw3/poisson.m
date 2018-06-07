function [uh, in] = poisson(f, fd, h0, p, t)
    geps = .001 * h0; 
    ind = (feval(fd,p) < -geps); % find interior nodes
    Np = size(p, 1); 
    N = sum(ind); % Np nodes; N is interior nodes
    in = zeros(Np, 1); 
    in(ind) = (1:N)'; % number the interior nodes
    ff = zeros(Np, 1);
    for j = 1:Np 
        ff(j) = feval(f, p(j, :)); 
    end % eval f once for each node
    % loop over triangles to set up stiffness matrix A and load vector b
    A = sparse(N, N); 
    b = zeros(N, 1);
    for n = 1:size(t,1)
        j = t(n, 1); 
        k = t(n, 2); 
        l = t(n, 3); 
        vj = in(j); 
        vk = in(k); 
        vl = in(l);
        J = [p(k,1)-p(j,1), p(l,1)-p(j,1); p(k,2)-p(j,2), p(l,2)-p(j,2)];
        ar = abs(det(J))/2; 
        C = ar/12; 
        Q = inv(J'*J); 
        fT = [ff(j) ff(k) ff(l)];
        if vj > 0
            A(vj,vj) = A(vj,vj) + ar*sum(sum(Q)); 
            b(vj) = b(vj) + C*fT*[2 1 1]'; 
        end
        if vk > 0
            A(vk,vk) = A(vk,vk) + ar*Q(1,1); 
            b(vk) = b(vk) + C*fT*[1 2 1]'; 
        end
        if vl > 0
            A(vl,vl) = A(vl,vl) + ar*Q(2,2); 
            b(vl) = b(vl) + C*fT*[1 1 2]'; 
        end
        if vj*vk > 0
            A(vj,vk) = A(vj,vk) - ar*sum(Q(:,1)); 
            A(vk,vj) = A(vj,vk); 
        end
        if vj*vl > 0
            A(vj,vl) = A(vj,vl) - ar*sum(Q(:,2)); 
            A(vl,vj) = A(vj,vl); 
        end
        if vk*vl > 0
           A(vk,vl) = A(vk,vl) + ar*Q(1,2); 
           A(vl,vk) = A(vk,vl);  
        end
    end
    uh = zeros(Np,1); 
    uh(ind) = A\b;
    trimesh(t,p(:,1),p(:,2),uh), axis tight
end