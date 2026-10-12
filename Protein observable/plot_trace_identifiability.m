function [] = plot_trace_identifiability()

%%% command:
%%% plot_trace()

clc;
close all
set(0,'DefaultFigureVisible','on');

if ~exist('./figures', 'dir')
    mkdir('./figures')
end

a{1} = '$\mathrm{Log}_{10}$ $k_\mathrm{1}$';
a{2} = '$\mathrm{Log}_{10}$ $k_\mathrm{2}$';
a{3} = '$\mathrm{Log}_{10}$ $k_\mathrm{3}$';
a{4} = '$\mathrm{Log}_{10}$ $k_\mathrm{4}$';

labs = {a{1},...
    a{2},...
    a{3},...
    a{4}};

col = {[0 1 0]; [0 0 1]; [1 0 1]; [1 1 0]};

S1_struct = load('../chains/chains.mat');      % first chain file
S2_struct = load('./chains/chains.mat');     % second chain file
C1 = S1_struct.chains; % numeric matrix iterations x parameters
C2 = S2_struct.chains;

tiledplot = tiledlayout(1,4,'TileSpacing','Compact');
set(gcf, 'Position',  [500, 100, 1200, 300]);
mm = 0;
for aa = 1:length(labs)
    mm = mm+1;
    ax(mm) = nexttile(mm);
    % Ensure both C1 and C2 are numeric and have at least aa columns
    if aa > size(C1,2) || aa > size(C2,2)
        error('Parameter index exceeds number of columns in C1 or C2.');
    end

    hold(ax(mm),'on');                       % enable overlaying plots
    it1 = 1:size(C1,1);
    it2 = 1:size(C2,1);

    y1 = log10(C1(:,aa));
    y2 = log10(C2(:,aa));

    % Plot chain 1 (mRNA and protein observable)
    plot(ax(mm), it1, y1, ...
        'Color', [0 0 0], 'LineWidth', 0.8, 'LineStyle', '-', ...
        'DisplayName', 'mRNA and protein observable');

    % Plot chain 2 (protein observable)
    plot(ax(mm), it2, y2, ...
        'Color', col{aa}, 'LineWidth', 0.8, 'LineStyle', '-', ...
        'DisplayName', 'protein observable');

    yAll = [y1; y2];
    yAll = yAll(isfinite(yAll));
    if~isempty(yAll)
        yMin = min(yAll);
        yMax = max(yAll);
        yPad = 0.40 * max(1, yMax-yMin);
        % Set y-axis limits based on the calculated min and max values
        ylim(ax(mm), [yMin, yMax + yPad]);
    end    

    legend(ax(mm), 'show', 'Location', 'north');
    ylabel(ax(mm), labs{aa}, 'FontSize', 12, 'Interpreter', 'latex');
    xlabel(ax(mm), 'Iteration', 'FontSize', 12);

    % Use a single xlim that covers both chains
    set(ax(mm), 'xlim', [1 max(length(it1), length(it2))]);

    hold(ax(mm),'off');

end

tiledplot.TileSpacing = 'compact';
tiledplot.Padding = 'compact';

savefig(strcat('figures/trace_ide.fig'));
exportgraphics(gcf,strcat('figures/trace_ide.png'));

end