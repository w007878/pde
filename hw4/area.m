function S = area(x1, y1, x2, y2, x3, y3)
    xA = x2 - x1; 
    yA = y2 - y1;
    xB = x3 - x1;
    yB = y3 - y1;
    S = 0.5 * abs(xA*yB-xB*yA);
end