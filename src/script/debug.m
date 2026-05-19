fname = argv(){1};
plname = argv(){2};
sdfname = argv(){3};

data = dlmread(fname, ',');

X    = data(:, 1:5:end);
Y    = data(:, 2:5:end);
Zi    = data(:, 3:5:end);
GX   = data(:, 4:5:end);
GY   = data(:, 5:5:end);


% log scale
Z = sign(Zi) .* log1p(abs(Zi));

minx = X(1,   1);
maxx = X(1, end);
miny = Y(1,   1);
maxy = Y(end, 1);

x = linspace(minx, maxx, size(Z, 2));
y = linspace(miny, maxy, size(Z, 1));


zmax = max(Z(:));
zmin = min(Z(:));

zfloor = zmax + 0.1;

negative_lev = linspace(zmin, 0, 10)(1:end-1);
positive_lev = linspace(0, zmax, 10)(2:end);
levels = [negative_lev, positive_lev];

fig = figure();

% SURFACE
surf(x, y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.5);
hold on;

% CONTOUR LINES
contour3(x, y, Z, levels, 'LineWidth', 1);
contour3(x, y, Z, [0 0], 'r', 'LineWidth', 2);

% GRADIENT FIELD
n = 6;
hGrad = quiver3(X(1:n:end, 1:n:end), Y(1:n:end, 1:n:end), zeros(size(X(1:n:end, 1:n:end))) + zfloor, GX(1:n:end, 1:n:end), GY(1:n:end, 1:n:end), zeros(size(GX(1:n:end, 1:n:end))), 0.8, 'b');

% POLYLINE
pl = dlmread(plname, ',');
hPoly = [];
for i = 1:size(pl, 1)
    h = plot3([pl(i,1), pl(i,3)], [pl(i,2), pl(i,4)], [zfloor, zfloor], 'g-', 'LineWidth', 2);
    hPoly = [hPoly, h];
end

% POLYLINE NORMALS
hNormals = [];
for i = 1:size(pl, 1)
    % edge midpoint
    mx = (pl(i,1) + pl(i,3)) / 2;
    my = (pl(i,2) + pl(i,4)) / 2;
    % edge direction
    dx = pl(i,3) - pl(i,1);
    dy = pl(i,4) - pl(i,2);

    nx = -dy;
    ny =  dx;
    len = sqrt(nx^2 + ny^2);
    nx = nx / len;
    ny = ny / len;

    h = quiver3(mx, my, zfloor, nx, ny, 0, 0.2, 'r', 'LineWidth', 1.5);
    hNormals = [hNormals, h];
end

% POINTS AND BETAS
sdfpts = dlmread(sdfname, ',');
npts = size(sdfpts, 1);
hSdfPts  = plot3(sdfpts(:,1), sdfpts(:,2), zeros(npts,1) + zfloor, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
hSdfVecs = quiver3(sdfpts(:,1), sdfpts(:,2), zeros(npts,1) + zfloor, sdfpts(:,3), sdfpts(:,4), zeros(npts,1), 0.3, 'k', 'LineWidth', 2);

colormap(cool);
colorbar;
axis equal;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');

view(0, 90);




