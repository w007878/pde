function phi = getPhi(px, py, i, x, y, nTri, tri, assTri)
    for j = 1:nTri
        if(assTri(i, j) && checkIn(px, py, x(tri(j, 1)), y(tri(j, 1)), x(tri(j, 2)), y(tri(j, 2)), x(tri(j, 3)), y(tri(j, 3))))

            if tri(j, 1) == i
                pB = tri(j, 2);
                pC = tri(j, 3);
            elseif tri(j, 2) == i
                pB = tri(j, 1);
                pC = tri(j, 3);
            else
                pB = tri(j, 2);
                pC = tri(j, 1);
            end

            para = [1; 0; 0] \ [x(i), y(i), 1; x(pB), y(pB), 1; x(pC), y(pC), 1];
            phi = para * [px; py; 1];
            return;
        end
    end
    phi = 0;
end