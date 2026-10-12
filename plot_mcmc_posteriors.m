function [] = plot_mcmc_posteriors()

%%% command:
%%% plot_mcmc_posteriors()

clc;
close all
set(0,'DefaultFigureVisible','on');

if ~exist('./figures', 'dir')
    mkdir('./figures')
end

chainStruct = load('./chains/chains.mat');

chains = chainStruct.chains(2:2:end,:);


%%% Aggregated posterior plots %%%
fontsize = 12;
leg = {'k_{1}',...
    'k_{2}',...
    'k_{3}',...
    'k_{4}'};

true_values = [20, 0.5, 0.25, 0.5];
%truecolor = {[0.6 0 0.6],[0 0.6 0],[0.6 0.6 0],[0 0.6 0.6]};

fcolor = {[0 0 0],[1 0 0],[0 0 1], [0 1 0]};
figure(1);
tiledplot = tiledlayout(1,4,'TileSpacing','Compact');
set(gcf, 'Position',  [300, 100, 800, 250]);
for aa = 1:4
    ax(aa) = nexttile(aa);            
    set(ax(aa),...
    'box','on',...
    'FontSize',10);
set(gca,'TickLength',[0.025, 0.01])
hold on;
histogram(log10(chains(:,aa)),...
         'FaceColor', 'none',...
         'EdgeColor', fcolor{aa},...
         'FaceAlpha', 1,...
         'EdgeAlpha', 1,...
         'Normalization', 'probability',...
         'DisplayStyle', 'stairs',...
         'LineWidth',1.5);
xlabel(ax(aa),strcat('Log_{10} ', {' '}, leg(aa)),'FontSize',fontsize);
yl = ylim();
xval = log10(true_values(aa));
 %plot vertical line spanning the full y range
 plot([xval xval], yl, 'Color', fcolor{aa}, 'LineWidth', 2, 'LineStyle', '--');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nn = 1:4
%     histogram(log10(chains(:,nn)),...
%         'FaceColor', 'none',...
%         'EdgeColor', fcolor{nn},...
%         'FaceAlpha', 1,...
%         'EdgeAlpha', 1,...
%         'Normalization', 'probability',...
%         'DisplayStyle', 'stairs',...
%         'LineWidth',1.5);
%     hold on;
% end
end
% % After histograms are drawn, get current y-limits to span the vertical lines
% yl = ylim();
% 
% % Plot true-value vertical lines
% for nn = 1:4
%     xval = log10(true_values(nn));
%     % plot vertical line spanning the full y range
%     plot([xval xval], yl, 'Color', truecolor{nn}, 'LineWidth', 2, 'LineStyle', '--');
% end

% l=legend(leg);
% set(l,'Location','northwest');
% set(l,'FontSize',fontsize);
% set(l,'interpreter','latex');

ylabel(tiledplot,'Normalized frequency','FontSize',fontsize);
tiledplot.TileSpacing = 'compact';
tiledplot.Padding = 'compact';

savefig(strcat('figures/posteriors.fig'));
exportgraphics(gcf,'figures/posteriors.png');
exportgraphics(gcf,'figures/posteriors.eps','ContentType','vector');               

end