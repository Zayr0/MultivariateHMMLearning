function plotCovarianceMatrices(O)
n = length(O);
figure('Renderer', 'painters', 'Position', [10 10 1500 1000])

for i = 1:n
    covMat = corrcoef(O{i});

    subplot(1, n, i);
    h = heatmap(covMat, Colormap=jet());
    h.Title = ['Correlation_' num2str(i) ' - Sum: ' num2str(sum(covMat, "all")/size(covMat, 1)^2)];
end
end