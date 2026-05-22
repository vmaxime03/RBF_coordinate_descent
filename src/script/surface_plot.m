folder = argv(){1};
fname        = fullfile(folder, 'sdf.csv');
plname       = fullfile(folder, 'polyline.csv');
sdfname      = fullfile(folder, 'sdf_params.csv');


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

negative_lev = linspace(zmin, 0, 10)(1:end-1);
positive_lev = linspace(0, zmax, 10)(2:end);
levels = [negative_lev, positive_lev];

fig = figure();

% SURFACE
surf(x, y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.5);
hold on;

% CONTOUR LINES
contour3(x, y, Z, levels, 'LineWidth', 1);
contour3(x, y, Z, [0 0], 'r', 'LineWidth', 4);

% POLYLINE
pl = dlmread(plname, ',');
for i = 1:size(pl, 1)
    h = plot3([pl(i,1), pl(i,3)], [pl(i,2), pl(i,4)], [0, 0], 'g-', 'LineWidth', 2, 'Clipping', 'off');
end
set(gca, 'SortMethod', 'childorder');

% CONTROL POINTS
sdfpts = dlmread(sdfname, ',');
npts = size(sdfpts, 1);
pts = plot3(sdfpts(:,1), sdfpts(:,2), zeros(npts,1), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'magenta');



colormap(cool);
colorbar;
axis equal;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');

view(0, 90);

pause()
