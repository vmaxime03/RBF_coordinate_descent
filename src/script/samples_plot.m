folder   = argv(){1};
plname   = fullfile(folder, 'polyline.csv');
errname  = fullfile(folder, 'samples_error.csv');

pl   = dlmread(plname,    ',');
errd = dlmread(errname,   ',');

% DATA
nsamples = size(errd, 1);

dist_vals  = errd(:, 3);
grad_vals  = errd(:, 4);

% Secure the ranges against microscopic floating-point differences
dist_range = max(dist_vals) - min(dist_vals);
if dist_range < 1e-9; dist_range = 1; end
dist_norm = (dist_vals - min(dist_vals)) / dist_range; % Normalized strictly 0 to 1
 
grad_range = max(grad_vals) - min(grad_vals);
if grad_range < 1e-9; grad_range = 1; end
grad_norm = (grad_vals - min(grad_vals)) / grad_range; % Normalized strictly 0 to 1
 

cmap_dist = hot(256);
cmap_grad = cool(256);

figure();
hold on;

% DISTANCE ERROR — normalized bars toward +z, hot colormap
for i = 1 : nsamples
  % dist_norm is already 0-1, so just scale to 255 here
  ci  = max(1, round(dist_norm(i) * 255) + 1);
  col = cmap_dist(ci, :);
  plot3([errd(i,1), errd(i,1)], [errd(i,2), errd(i,2)], [0, dist_norm(i)], ...
        '-', 'Color', col, 'LineWidth', 1.5);
end
 
% GRADIENT ERROR — normalized bars toward -z, cool colormap
for i = 1 : nsamples
  % grad_norm is already 0-1, so just scale to 255 here
  ci  = max(1, round(grad_norm(i) * 255) + 1);
  col = cmap_grad(ci, :);
  plot3([errd(i,1), errd(i,1)], [errd(i,2), errd(i,2)], [0, -grad_norm(i)], ...
        '-', 'Color', col, 'LineWidth', 1.5);
end
 
% POLYLINE
for i = 1 : size(pl, 1)
  plot3([pl(i,1), pl(i,3)], [pl(i,2), pl(i,4)], [0, 0], 'g-', 'LineWidth', 2);
end
 

h = legend({sprintf("Dist err [%.4f, %.4f]", min(dist_vals), max(dist_vals)), ... 
            sprintf("Grad err [%.4f, %.4f]", min(grad_vals), max(grad_vals))}, ... 
            "location", "east"
    );
set(h, "fontsize", 20);

% Note: If bars still look completely flat, comment out 'axis equal' 
% to let the Z-axis scale independently from the X/Y plane.
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z (Normalized Error)');
view(30, 30);
 
pause();
