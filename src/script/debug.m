fname = argv(){1};
plname = argv(){2};

data = dlmread(fname, ',');

w = size(data, 1);
h = size(data, 2) / 5;

X    = data(:, 1:5:end);
Y    = data(:, 2:5:end);
% distance
Z    = data(:, 3:5:end);

% gradient
GX   = data(:, 4:5:end);
GY   = data(:, 5:5:end);


minx = X(1,   1);
maxx = X(1, end);
miny = Y(1,   1);
maxy = Y(end, 1);

x = linspace(minx, maxx, size(Z, 1));
y = linspace(miny, maxy, size(Z, 2));

Zt = sign(Z) .* log1p(abs(Z)); % pour que les ligne de niveaux ne soit pas trop eloignées

zmax = max(Zt(:));
zmin = min(Zt(:));

negative_lev = linspace(zmin, 0, 10)(1:end-1);
positive_lev = linspace(0, zmax, 10)(2:end);

levels =  [negative_lev, positive_lev];


% SURFACE PLOT
fig1 = figure();

surf(Z);

xlabel("x");
ylabel("y");
zlabel("f(x, y)");




% CONTOUR PLOT
fig2 = figure();

contour(x, y, Zt, levels);

hold on;
contour(x, y, Zt, [0 0], "r", "LineWidth", 2);

colorbar;
axis equal;

colormap(cool);




% GRADIENT FIELD

norm_G = sqrt(GX.^2 + GY.^2);
norm_G_safe = max(norm_G, 1e-12);

GXn = GX ./ norm_G_safe;
GYn = GY ./ norm_G_safe;

% quiver(X, Y, GXn, GYn, 0.8);

n = 6;
quiver(X(1:n:end, 1:n:end), Y(1:n:end, 1:n:end), GX(1:n:end, 1:n:end), GY(1:n:end, 1:n:end), 0.8);


hold on;
contour(x, y, Zt, [0 0], "r", "LineWidth", 6);



pl = dlmread(plname, ',');
for i = 1:size(pl, 1)
  plot([pl(i,1), pl(i,3)], [pl(i,2), pl(i,4)], 'g-', 'LineWidth', 2);
end



axis equal;
xlabel("x"); ylabel("y");

pause();
