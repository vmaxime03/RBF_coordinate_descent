folder       = argv(){1};
plname       = fullfile(folder, 'polyline.csv');
sdfname      = fullfile(folder, 'sdf_params.csv');
samplesnames = fullfile(folder, 'samples.csv');
errname      = fullfile(folder, 'samples_error.csv');

pl      = dlmread(plname,       ',');
sdfpts  = dlmread(sdfname,      ',');
samples = dlmread(samplesnames, ',');
errd    = dlmread(errname,      ',');

% DATA
nsamples = size(errd, 1);
npts     = size(sdfpts, 1);

dist_vals  = errd(:, 3);
grad_vals  = errd(:, 4);


dist_range = max(dist_vals) - min(dist_vals);
if dist_range == 0; dist_range = 1; end
dist_norm = (dist_vals - min(dist_vals)) ;% / dist_range;
 
grad_range = max(grad_vals) - min(grad_vals);
if grad_range == 0; grad_range = 1; end
grad_norm = (grad_vals - min(grad_vals)) ; % / grad_range;
 

cmap_dist = hot(256);
cmap_grad = cool(256);

figure();
hold on;

% DISTANCE ERROR — normalised bars toward +z, hot colormap
for i = 1 : nsamples
  ci  = max(1, round(dist_norm(i) * 255) + 1);
  col = cmap_dist(ci, :);
  plot3([errd(i,1), errd(i,1)], [errd(i,2), errd(i,2)], [0, dist_norm(i)], ...
        '-', 'Color', col, 'LineWidth', 1.5);
end
 
% GRADIENT ERROR — normalised bars toward -z, cool colormap
for i = 1 : nsamples
  ci  = max(1, round(grad_norm(i) * 255) + 1);
  col = cmap_grad(ci, :);
  plot3([errd(i,1), errd(i,1)], [errd(i,2), errd(i,2)], [0, -grad_norm(i)], ...
        '-', 'Color', col, 'LineWidth', 1.5);
end
 
% POLYLINE
for i = 1 : size(pl, 1)
  plot3([pl(i,1), pl(i,3)], [pl(i,2), pl(i,4)], [0, 0], 'g-', 'LineWidth', 2);
end
 
% SDF CONTROL POINTS + BETA VECTORS
plot3(sdfpts(:,1), sdfpts(:,2), zeros(npts,1), ...
      'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
quiver3(sdfpts(:,1), sdfpts(:,2), zeros(npts,1), ...
        sdfpts(:,3), sdfpts(:,4), zeros(npts,1), ...
        0.3, 'k', 'LineWidth', 2);

h = legend({sprintf("Dist err [%.4f, %.4f]", min(dist_vals), max(dist_vals)), ... 
			sprintf("Grad err [%.4f, %.4f]", min(grad_vals), max(grad_vals))}, ... 
			"location", "east"
	);
set(h, "fontsize", 20);

axis equal; grid on;
xlabel('x'); ylabel('y');
view(30, 30);
 
pause();
 

