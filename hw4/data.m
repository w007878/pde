% This script generate an ellipse and n nodes on the boundary, m nodes inside it.
n = 20;
m = 50;
eps = 1e-8;

%% The first part: generate an ellipse.
a = 2;
e = 0.5;
c = a * e;
b = sqrt(a^2 - c^2);

theta = 2*pi/n .* (1:n);
rho = b ./ sqrt(1 - e^2*power(cos(theta),2));

xBorder = rho .* cos(theta);
yBorder = rho .* sin(theta);

%% The second part: generate the inside points.
thetaIns = random('Uniform', 0, 2*pi, m, 1)';
rhoIns = random('Uniform', 0.3+eps, b./sqrt(1 - e^2*power(cos(thetaIns),2))-eps);
xIns = rhoIns .* cos(thetaIns);
yIns = rhoIns .* sin(thetaIns);

plot(xBorder, yBorder, '-', xIns, yIns, '.')

x = [xBorder xIns];
y = [yBorder yIns];

%% The third part: Write the points into a file and call triangle program
polyFile = fopen('points.poly', 'w');
fprintf(polyFile, '%d 2 1 0\n', n+m);
for i = 1:n
    fprintf(polyFile, '%d %.2f %.2f 0\n', i, xBorder(i), yBorder(i));
end
for i = 1:m
    fprintf(polyFile, '%d %.2f %.2f 0\n', i+n, xIns(i), yIns(i));
end
fprintf(polyFile, '%d 0\n', n);
for i = 1:n
    fprintf(polyFile, '%d %d %d\n', i, i, mod(i,n)+1);
end
fprintf(polyFile, '0\n');
fclose(polyFile);

%% Read file from the result of triangle
tri = fopen('points.1.ele');
tmp = fscanf(tri, '%d');
nTri = tmp(1);
triangle = reshape(tmp(4:size(tmp)),4, nTri)';
triangle = triangle(:, 2:4);
ass = zeros(n+m, n+m);
xCent = zeros(nTri);
yCent = zeros(nTri);

hold on
for i = 1:nTri
    xtmp = x(triangle(i,:));
    ytmp = y(triangle(i,:));
    plot([xtmp xtmp(1)], [ytmp ytmp(1)])
    ass(triangle(i,1), triangle(i,2)) = 1;
    ass(triangle(i,1), triangle(i,3)) = 1;
    ass(triangle(i,2), triangle(i,1)) = 1;
    ass(triangle(i,2), triangle(i,3)) = 1;
    ass(triangle(i,3), triangle(i,1)) = 1;
    ass(triangle(i,3), triangle(i,2)) = 1;    
    xCent(i) = mean(x(triangle(i, :)));
    yCent(i) = mean(y(triangle(i, :)));
end

assTri = zeros(n + m, nTri);
for i=1:nTri
    assTri(triangle(i, 1), i) = 1;
    assTri(triangle(i, 2), i) = 1;
    assTri(triangle(i, 3), i) = 1;
end