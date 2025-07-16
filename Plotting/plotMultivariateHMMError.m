function plotMultivariateHMMError(T, O, pi, T_hat, O_hats, pi_hat)
    if(length(O) ~= length(O_hats))
        warning("O and O_hats are not the same length");
        return;
    end

    n = min(length(O), 5);

    figure('Renderer', 'painters', 'Position', [100 100 2000 1500])

    subplot(3, n+2, 1);
    h = heatmap(pi, Colormap=jet());
    h.Title = '\pi';

    subplot(3, n+2, 2);
    h = heatmap(T, Colormap=jet());
    h.Title = 'T';
    
    for i = 1:n
        subplot(3, n+2, i+2);
        h = heatmap(O{i}, Colormap=jet());
        h.Title = ['O_' num2str(i)];
    end

    subplot(3, n+2, n+3);
    h = heatmap(pi_hat, Colormap=jet());
    h.Title = '\pi^{hat}';


    subplot(3, n+2, n+4);
    h = heatmap(T_hat, Colormap=jet());
    h.Title = 'T^{hat}';
    
    for i = 1:n
        subplot(3, n+2, (n+3)+i+1);
        h = heatmap(O_hats{i}, Colormap=jet());
        h.Title = ['O^{hat}_' num2str(i)];
    end

    subplot(3, n+2, 2*n+5);
    h = heatmap(abs(pi - pi_hat), Colormap=jet());
    h.Title = '\pi^{error}';

    subplot(3, n+2, 2*n+6);
    h = heatmap(abs(T - T_hat), Colormap=jet());
    h.Title = 'T^{Error}';
    
    for i = 1:n
        subplot(3, n+2, 2*(n+2)+i+2);
        h = heatmap(abs(O{i} - O_hats{i}), Colormap=jet());
        h.Title = ['O^{Error}_' num2str(i)];
    end
    
end