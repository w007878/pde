function [p, t] = recmesh2d(gap, con)
    nx = (con(2, 1)-con(1, 1))/gap + 1;
    ny = (con(2, 2)-con(1, 2))/gap + 1;
    disp([nx, ny]);
    n = nx * ny;
    p = zeros(nx * ny, 2);
    t = [];
    label = @(i, j)(i-1)*nx+j;
    
    for i = 1:nx
        for j = 1:ny
            p(label(i,j),:) = [gap*(i-1)+con(1,1), gap*(j-1)+con(1,2)];
            if i > 1 && j < ny
                t = [t; label(i,j), label(i-1,j), label(i,j+1)];
            end
        end
    end
    disp(p);
    disp(t);
end