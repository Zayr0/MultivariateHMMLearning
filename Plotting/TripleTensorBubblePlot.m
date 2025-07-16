function TripleTensorBubblePlot(v1, v2, d)

% Sample 3D data
x = kron(ones(1, d), kron(ones(1, d), 1:d));
y = kron(ones(1, d), kron(1:d, ones(1, d)));
z = kron(kron(1:d, ones(1, d)), ones(1, d));

v_diff = v1 - v2;       % Difference between tensors

% Set the colormap and ensure a common color scale
cmap = parula; % Choose a colormap
clim = [min([v1; v2; v_diff]), max([v1; v2; v_diff])]; % Common color limits


% Create a figure
figure('Renderer', 'painters', 'Position', [10 10 1500 500]);

% First subplot (Tensor 1)
subplot(1,3,1);
bubblechart3(x, y, z, abs(v1), v1); % Bubble size and color based on v1
caxis(clim); % Set common color limits
colormap(cmap); % Apply colormap
title('P3');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; hold on;

% Second subplot (Tensor 2)
subplot(1,3,2);
bubblechart3(x, y, z, abs(v2), v2); % Bubble size and color based on v2
caxis(clim); % Maintain color limits
colormap(cmap);
title('M');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; hold on;

% Third subplot (Difference)
subplot(1,3,3);
bubblechart3(x, y, z, abs(v_diff), abs(v_diff)); % Bubble size and color based on difference
caxis(clim);
colormap(cmap);
title('Difference');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; hold on;

% Create a shared colorbar
cb = colorbar;
% cb.Layout.Tile = 'east'; % Place colorbar on the right side
% cb.Label.String = 'Value (Color Scale)';

% Apply the same colormap across all subplots
colormap(cmap);

end