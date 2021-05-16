fprintf('loading, plotting and save opt history ... \n')
load('test_history')
fig = figure('visible','off');
plot_opt_history(G_history.iter,G_history.fval,G_history.cstr(:,1))
figName = ['test_history'];
set(fig, 'PaperPositionMode', 'auto');
print('-djpeg ', '-r300', figName);