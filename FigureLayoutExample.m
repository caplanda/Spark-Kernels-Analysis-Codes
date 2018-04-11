clear; clc; close all;

figure( 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6], ...
        'Color', 'white' ) ;
    % - Bg axes and main title.
    bgAxes = axes( 'Position', [0, 0, 1, 1], 'XColor', 'none', 'YColor', 'none', ...
        'XLim', [0, 1], 'YLim', [0, 1] ) ;
    text( 0.5, 0.95, 'A Small Example', 'FontSize', 16, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold' ) ;
    % - Positions.
    x = linspace( 0.05, 0.6, 4 ) ;
    w = 0.8 * diff( x(1:2) ) ;
    y = linspace( 0.07, 0.7, 3 ) ;
    h = 0.8 * diff( y(1:2) ) ;
    % - Headers for array still in bgAxes.
    for colId = 1 : numel( x )-1
        text( x(colId)+w/2, 0.03, sprintf( 'Label X%d', colId ), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold' ) ;
    end
    for rowId = 1 : numel( y )-1
        text( 0.02, y(rowId)+h/2, sprintf( 'Label Y%d', rowId ), ...
            'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold' ) ;
    end
    % - Array of plots.
    for colId = 1 : numel( x )-1
        for rowId = 1 : numel( y )-1
            axes( 'Position', [x(colId), y(rowId), w, h] ) ;
            plot( sin( rand(1) * (1:10)), 'b' ) ;
            grid( 'on' ) ;
        end
    end
    % - Surface.
    axes( 'Position', [0.65, 0.45, 0.3, 0.45] ) ;
    [X, Y] = meshgrid( -5: .5 : 5 ) ;
    Z = Y.*sin(X) - X.*cos(Y) ;
    s = surf(X,Y,Z,'FaceAlpha',0.5) ;
    s.EdgeColor = 'none';
    % - Time series
    axes( 'Position', [0.05, 0.72, 0.515, 0.15] ) ;
    x = linspace( 0, 10, 100 ) ;
    plot( rand(size(x)) +  3 * sin(x) .* exp(-x/5), 'r' ) ;
    xlabel( 't [s]' ) ;
    ylabel( 'A [V]' ) ;
    grid( 'on' ) ;
    % - Barchart.
    axes( 'Position', [0.65, 0.05, 0.3, 0.3] ) ;
    histogram( randn(1000, 1) ) ;
    set( gca, 'Box', 'off' ) ;