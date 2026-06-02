% clusters_plot.m
% Usage: octave clusters_plot.m <output_folder> <iteration>
% Example: octave clusters_plot.m build/test/ 3

folder    = argv(){1};
iteration = str2num(argv(){2});

plname   = fullfile(folder, 'polyline.csv');
sdfname  = fullfile(folder, 'sdf_params.csv');
clname   = fullfile(folder, sprintf('%d_clusters.csv', iteration));

pl     = dlmread(plname,  ',');
sdfpts = dlmread(sdfname, ',');

% cluster file: cluster_id, x, y, dist_err, grad_err
cl = [];
if exist(clname, 'file')
  cl = dlmread(clname, ',');
end

% ---- colour palette: one colour per cluster --------------------------------
nclust = 0;
if ~isempty(cl)
  nclust = max(cl(:, 1)) + 1;   % cluster ids are 0-based
end

cmap = hsv(max(nclust, 1));

% ---- normalise errors (keep absolute scale for legend) --------------------
dist_vals = [];
grad_vals = [];
if ~isempty(cl)
  dist_vals = cl(:, 4);
  grad_vals = cl(:, 5);
end

dist_max = max(abs(dist_vals));
grad_max = max(abs(grad_vals));
if isempty(dist_max) || dist_max == 0; dist_max = 1; end
if isempty(grad_max) || grad_max == 0; grad_max = 1; end

dist_norm = dist_vals / dist_max;
grad_norm = grad_vals / grad_max;

figure();
hold on;

% ---- polyline --------------------------------------------------------------
for i = 1 : size(pl, 1)
  plot3([pl(i,1), pl(i,3)], [pl(i,2), pl(i,4)], [0, 0], ...
        'g-', 'LineWidth', 2);
end

% ---- cluster samples — dist error +z, grad error -z, coloured by cluster --
for i = 1 : size(cl, 1)
  cid = cl(i, 1) + 1;          % 1-based
  col = cmap(cid, :);
  x   = cl(i, 2);
  y   = cl(i, 3);

  % distance error bar upward
  plot3([x, x], [y, y], [0,  dist_norm(i)], ...
        '-',  'Color', col, 'LineWidth', 2);

  % gradient error bar downward
  plot3([x, x], [y, y], [0, -grad_norm(i)], ...
        '--', 'Color', col, 'LineWidth', 2);
end

% ---- SDF control points + beta vectors -------------------------------------
npts = size(sdfpts, 1);
plot3(sdfpts(:,1), sdfpts(:,2), zeros(npts,1), ...
      'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
quiver3(sdfpts(:,1), sdfpts(:,2), zeros(npts,1), ...
        sdfpts(:,3), sdfpts(:,4), zeros(npts,1), ...
        0.3, 'k', 'LineWidth', 2);

% ---- legend ----------------------------------------------------------------
leg = {};
for ci = 0 : nclust-1
  leg{end+1} = sprintf('Cluster %d  (dist +z solid, grad -z dashed)', ci);
end

h = legend(leg, 'Location', 'eastoutside', 'FontSize', 12);

title(sprintf('Clusters — iter %d   dist [%.4f, %.4f]   grad [%.4f, %.4f]', ...
      iteration, min(dist_vals), max(dist_vals), min(grad_vals), max(grad_vals)));
xlabel('x'); ylabel('y'); zlabel('normalised error');
axis equal; grid on;
view(30, 30);

pause();
