fname = argv(){1};

% data = readmatrix(fname); % matlab 
data = dlmread(fname, ',');

w = size(data, 1);
h = size(data, 2) / 3;

X = data(:, 1:3:end);
Y = data(:, 2:3:end);
Z = data(:, 3:3:end);

minx = X(1,   1);
maxx = X(end, 1);
miny = Y(1,   1);
maxy = Y(1,   end);

x = linspace(minx, maxx, size(Z, 1));
y = linspace(miny, maxy, size(Z, 2));

Zt = sign(Z) .* log1p(abs(Z)); % pour que les ligne de niveaux ne soit pas trop eloignées

zmax = max(Zt(:));
zmin = min(Zt(:));

negative_lev = linspace(zmin, 0, 10)(1:end-1);
positive_lev = linspace(0, zmax, 10)(2:end);

levels =  [negative_lev, positive_lev];


fig1 = figure();

surf(Z);

% hold on;
% contour3(x, y, Z, [0 0], "-r", 'LineWidth', 5);

xlabel("x");
ylabel("y");
zlabel("f(x, y)");



fig2 = figure();

contour(x, y, Zt, levels);

hold on;
contour(x, y, Zt, [0 0], "r", "LineWidth", 2);

colorbar;
axis equal;

colormap(cool);

xlabel("x");
ylabel("y");

pause();
