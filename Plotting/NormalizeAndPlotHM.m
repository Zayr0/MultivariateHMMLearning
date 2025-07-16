function NormalizeAndPlotHM(A, B, C)
i = 1;

if (nargin < 3)
    B = Stochasticize(A{2});
    C = Stochasticize(A{3});
    A = Stochasticize(A{1});
else
    A = Stochasticize(A);
    B = Stochasticize(B);
    C = Stochasticize(C);
end

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

subplot(1, 3, 1);
h = heatmap(A, Colormap=jet());
h.Title = 'A';

subplot(1, 3, 2);
h = heatmap(B, Colormap=jet());
h.Title = 'B';

subplot(1, 3, 3);
h = heatmap(C, Colormap=jet());
h.Title = 'C';
end