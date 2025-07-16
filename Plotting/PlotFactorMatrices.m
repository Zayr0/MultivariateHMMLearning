function PlotFactorMatrices(U)
R = length(U);
figure('Renderer', 'painters', 'Position', [10 10 1200 600])

for i = 1:R
    subplot(1, R, i);
    h = heatmap(U{i}, Colormap=jet(20), FontSize=10);
    h.Title = ['B_' num2str(i)];
end
end