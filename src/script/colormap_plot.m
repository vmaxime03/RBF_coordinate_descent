folder   = argv(){1};
fname    = fullfile(folder, 'sdf.csv');
plname   = fullfile(folder, 'polyline.csv');
sdfname  = fullfile(folder, 'sdf_params.csv');

data = dlmread(fname, ',');

X  = data(:, 1:5:end);
Y  = data(:, 2:5:end);
Zi = data(:, 3:5:end);
GX = data(:, 4:5:end);
GY = data(:, 5:5:end);

% Z = sign(Zi) .* log1p(abs(Zi));
Z = Zi;

minx = X(1,   1);
maxx = X(1, end);
miny = Y(1,   1);
maxy = Y(end, 1);

x = linspace(minx, maxx, size(Z, 2));
y = linspace(miny, maxy, size(Z, 1));

fig = figure();

imagesc(x, y, Z); 
set(gca, 'YDir', 'normal'); 
colormap(cool);
colorbar;
hold on; 

contour(x, y, Z, [0 0], 'r', 'LineWidth', 2);

n = 6;
X_sub  = X(1:n:end, 1:n:end);
Y_sub  = Y(1:n:end, 1:n:end);
GX_sub = GX(1:n:end, 1:n:end);
GY_sub = GY(1:n:end, 1:n:end);

hGrad = quiver(X_sub, Y_sub, GX_sub, GY_sub, 0.8, 'k');

try
    pl = dlmread(plname, ',');
    for i = 1:size(pl, 1)
        plot([pl(i,1), pl(i,3)], [pl(i,2), pl(i,4)], 'g-', 'LineWidth', 2);
    end
catch
    disp('No polyline file found or failed to read.');
end

axis equal;
axis([minx maxx miny maxy]); 
xlabel('x'); 
ylabel('y');

hold off;
pause();
