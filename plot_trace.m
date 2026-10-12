function [] = plot_trace()

%%% command:
%%% plot_trace()

clc;
close all
set(0,'DefaultFigureVisible','on');

if ~exist('./figures', 'dir')
    mkdir('./figures')
end

a{1} = '$k_\mathrm{1}$';
a{2} = '$k_\mathrm{2}$';
a{3} = '$k_\mathrm{3}$';
a{4} = '$k_\mathrm{4}$';

labs = {a{1},...
    a{2},...
    a{3},...
    a{4}};

col = {[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [1 1 0]};

tiledplot = tiledlayout(1,4,'TileSpacing','Compact');
set(gcf, 'Position',  [500, 100, 1200, 300]);

    chains = load(strcat('./chains/chains.mat'));
    chains = chains.chains(:,:);
    mm = 0;
    for aa = 1:length(labs)
        mm = mm+1;
        ax(mm) = nexttile(mm);
        %plot(ax(mm),1:1:length(chains(:,aa)),log10(chains(:,aa)),'Color',[col 1],'LineWidth',0.1,'LineStyle','-');
        plot(ax(mm), 1:length(chains(:,aa)), log10(chains(:,aa)), 'Color', col{aa}, 'LineWidth', 0.1, 'LineStyle', '-');

        hold on;
        ylabel(ax(mm),labs{aa},'FontSize',12,'Interpreter','latex');
        xlabel(ax(mm),'Iteration','FontSize',12);
        set(ax(mm),'xlim',[1 length(chains(:,aa))]);
    end
  
tiledplot.TileSpacing = 'compact';
tiledplot.Padding = 'compact';

savefig(strcat('figures/trace.fig'));
exportgraphics(gcf,strcat('figures/trace.png'));
exportgraphics(gcf,strcat('figures/trace.eps'),'ContentType','vector');

end