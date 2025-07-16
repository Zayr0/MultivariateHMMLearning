function plotMultivariateMatricesAndSequences(seq, states, T, O, m)

n = length(O);
[d, k] = size(O{1});

seq = seq(:, 1:m);
states = states(1:m);

figure();
split = 1:10;
cmap = parula(k);  % state color map (no black, just state colors)
state_colors_rgb = cmap;  % k-by-3, each row is the color for a state
state_colors = states + 1; 
t = tiledlayout(n + 1, split(end), 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot observation subplots
for i = 1:n
    % Emission matrix subplot
    % subplot(n+2, split(end), (i-1) * split(end) + split(1));
    nexttile((i-1)*split(end) + 1);
    Oi = O{i};  % d x k
    
    % Brighten using sqrt scaling
    Oi_bright = Oi ./ max(Oi, [], 1);  % Apply brightness correction
    
    % Create RGB image
    [d, k] = size(Oi_bright);
    O_rgb = zeros(d, k, 3);
    for col = 1:k
        color = state_colors_rgb(col, :);
        for ch = 1:3
            O_rgb(:, col, ch) = Oi_bright(:, col) * color(ch);
        end
    end
    image(O_rgb);
    % title(["$O^{(", num2str(i), ")} = P(y_t|x_t)$"], "Interpreter", "latex", "FontSize", 15);
    ylabel("$D$", "Interpreter", "latex", "FontSize", 15);
    xlabel("$s_t$", "Interpreter", "latex", "FontSize", 15);
    set(gca, 'XTick', 1:k);

    % Sequence data colored by state
    % subplot(n+2, split(end), (i-1) * split(end) + split(2:end));
    nexttile((i-1)*split(end) + 2, [1 split(end)-1]);
    row_indices = seq(i, :);
    col_indices = 1:m;
    values = state_colors(:);
    A = full(sparse(row_indices, col_indices, values, d, m));
    s = pcolor(A);
    s.LineStyle = 'none';
    title(['Channel '  num2str(i)], "FontSize", 15);
    set(gca, 'YDir', 'reverse');
    colormap([0,0,0; cmap]);  % include black
    caxis([1 k+1]);
    grid on;
end

% --- Transition matrix with gradient coloring ---
% subplot(n+2, split(end), (n) * split(end) + split(1));
nexttile(n*split(end)+1);
T_bright = sqrt(sqrt(T));  % Brighten low values
T_rgb = zeros(k, k, 3);
for i = 1:k
    for j = 1:k
        color = state_colors_rgb(j, :);  % Use destination state's color
        for ch = 1:3
            T_rgb(i, j, ch) = T_bright(i, j) * color(ch);
        end
    end
end
image(T_rgb);
ylabel("$s_{t+1}$", "Interpreter", "latex", "FontSize", 15);
xlabel("$s_t$", "Interpreter", "latex", "FontSize", 15);
set(gca, 'XTick', 1:k, 'YTick', 1:k);

% --- State timeline plot ---
% subplot(n+2, split(end), (n) * split(end) + split(2:end));
nexttile(n*split(end)+2, [1 split(end)-1]);
imagesc(1:m, [0 1], state_colors);
colormap([0,0,0; cmap]);
clim([1 k+1]);
title("States", "FontSize", 15);
set(gca, 'YTickLabel', []);
xlabel("Timesteps", "Interpreter", "latex", "FontSize", 15);
hold on;



labels = arrayfun(@(s) ['State ' num2str(s)], 1:k, 'UniformOutput', false);
legend_handles = zeros(1,k);
hold on;
for s = 1:k
    legend_handles(s) = plot(NaN, NaN, 's', ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', cmap(s,:), ...
        'MarkerSize', 12);  % Increased marker size
end
hold off;

lgd = legend(legend_handles, labels, 'Location', 'southoutside', ...
    'Orientation', 'horizontal');
lgd.Box = 'off';
lgd.FontSize = 17;           % Increased font size
lgd.NumColumns = k;          % Optional: one column per state for better spacing
end

