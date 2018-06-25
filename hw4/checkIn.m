function res = checkIn(x, y, xA, yA, xB, yB, xC, yC)
    eps = 1e-8;
    area0 = area(xA, yA, xB, yB, xC, yC);
    area1 = area(xA, yA, xB, yB, x, y);
    area2 = area(xA, yA, x, y, xC, yC);
    area3 = area(x, y, xB, yB, xC, yC);
    res = abs(area0 - area1 - area2 - area3) < eps;
end
