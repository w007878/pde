function [z] = solve(n, m, x, y, nTri, tri, assTri, xCent, yCent, points, f)
    a = zeros(m, m);
    b = zeros(m, 1);
    
    for i = (n+1):(n+m)
        for j = 1:nTri
            if assTri(i, j)
                O = i;
                if tri(j, 1) == i 
                    A = tri(j, 2); B = tri(j, 3);
                elseif tri(j, 2) == i
                    A = tri(j, 1); B = tri(j, 3);
                else
                    A = tri(j, 1); B = tri(j, 2);
                end
                if (x(A)-x(O))*(y(B)-y(O)) - (y(A)-y(O))*(x(B)-x(O)) < 0
                    tmp = A; A = B; B = tmp;
                end
                if A <= n || B <= n 
                    continue; 
                end
                S = area(x(O), y(O), x(A), y(A), x(B), y(B));
                a(i-n, A-n) = a(i-n, A-n) + ((x(O)-x(B))*(x(B)-x(A)) + (y(O)-y(B))*(y(B)-y(A))) / (4*S);
                a(i-n, B-n) = a(i-n, B-n) + ((x(A)-x(O))*(x(B)-x(A)) + (y(A)-y(O))*(y(B)-y(A))) / (4*S);
                a(i-n, O-n) = a(i-n, O-n) + (-(x(O)-x(B))*(x(B)-x(A)) - (x(A)-x(O))*(x(B)-x(A)) -(y(O)-y(B))*(y(B)-y(A)) - (y(A)-y(O))*(y(B)-y(A))) / (4*S);
                               
                xC = xCent(j); yC = yCent(j);
                b(i-n) = b(i-n) + integral2(@(p,q)f(p*0.5*(x(A)-x(O))+q*(xC-x(O)), p*0.5*(y(A)-y(O))+q*(yC-y(O))), 0, 1, 0, @(x)1-x);
                b(i-n) = b(i-n) + integral2(@(p,q)f(p*0.5*(x(B)-x(O))+q*(xC-x(O)), p*0.5*(y(B)-y(O))+q*(yC-y(O))), 0, 1, 0, @(x)1-x);
            end
        end
    end
    
    u = b\a;
    %disp(u);
    
    z = zeros(size(points, 1), 1);
    phi = @(i, px, py)getPhi(px, py, i+n, x, y, nTri, tri, assTri);
    
    for i = 1:size(points, 1)
        for j = 1:m
            z(i) = z(i) + u(j) * phi(j, points(i, 1), points(i, 2));
%            fprintf('%f \n',u(j) *  phi(j, points(i), points(j)));
        end
    end

    plot3(points(:, 1), points(:, 2), z, '.');
end