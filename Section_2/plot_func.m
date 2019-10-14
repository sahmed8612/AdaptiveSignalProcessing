function plot_func(title_str, x, y, axis_labels, multiplier)
    axis(axis_labels)
    set(gca,'fontsize',20);
    title(title_str,'FontSize',20*multiplier);
    xlabel(x, 'FontSize', 20);
    ylabel(y, 'FontSize', 20);
end