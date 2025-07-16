function plotUnivariateHMMError(T, O, pi, T_hat, O_hat, pi_hat)

    figure('Renderer', 'painters', 'Position', [10 10 1500 900])
    
    if(nargin < 5)
    
        subplot(2, 3, 1);
        h = heatmap(T, Colormap=jet());
        h.Title = 'T';
        
        subplot(2, 3, 2);
        h = heatmap(T_hat, Colormap=jet());
        h.Title = 'T_{hat}';
        
        subplot(2, 3, 3);
        h = heatmap(abs(T-T_hat), Colormap=jet());
        h.Title = 'T_{error}';
        
        subplot(2, 3, 4);
        h = heatmap(O, Colormap=jet());
        h.Title = 'O';
        
        subplot(2, 3, 5);
        h = heatmap(O_hat, Colormap=jet());
        h.Title = 'O_{hat}';
        
        subplot(2, 3, 6);
        h = heatmap(abs(O-O_hat), Colormap=jet());
        h.Title = 'O_{error}';
    else
        subplot(3, 3, 1);
        h = heatmap(T, Colormap=jet());
        h.Title = 'T';
        
        subplot(3, 3, 2);
        h = heatmap(T_hat, Colormap=jet());
        h.Title = 'T_{hat}';
        
        subplot(3, 3, 3);
        h = heatmap(abs(T-T_hat), Colormap=jet());
        h.Title = 'T_{error}';
        
        subplot(3, 3, 4);
        h = heatmap(O, Colormap=jet());
        h.Title = 'O';
        
        subplot(3, 3, 5);
        h = heatmap(O_hat, Colormap=jet());
        h.Title = 'O_{hat}';
        
        subplot(3, 3, 6);
        h = heatmap(abs(O-O_hat), Colormap=jet());
        h.Title = 'O_{error}';
        
        subplot(3, 3, 7);
        h = heatmap(pi, Colormap=jet());
        h.Title = '\pi';
        
        subplot(3, 3, 8);
        h = heatmap(pi_hat, Colormap=jet());
        h.Title = '\pi_{hat}';
        
        subplot(3, 3, 9);
        h = heatmap(abs(pi-pi_hat), Colormap=jet());
        h.Title = '\pi_{error}';
    end
end