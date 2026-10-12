function [] = plot_practical_identifiability()

clc;
clear all;
close all;

if ~exist('./figures','dir')
    mkdir('./figures')
end

%% Load profile likelihood results
k_1 = load('./k_1.mat');
k_2 = load('./k_2.mat');
k_3 = load('./k_3.mat');
k_4 = load('./k_4.mat');

data = {k_1, k_2, k_3, k_4};

%% Load maximum likelihood values
MLvalue = load('./MLvalue.mat');
MLpars  = load('./MLpars.mat');

MLvalue = MLvalue.MLvalue;
MLpars  = MLpars.MLpars;

%% Labels
figlabs = {'(A)','(B)','(C)','(D)'};

xLabels = {...
    '$k_1\,(\mathrm{h}^{-1})$',...
    '$k_2\,(\mathrm{h}^{-1})$',...
    '$k_3\,(\mathrm{h}^{-1})$',...
    '$k_4\,(\mathrm{h}^{-1})$'};

fontsize = 14;

%% 95% confidence threshold
threshold = -MLvalue + 0.5*chi2inv(0.95,1);

%% Plot
tiledplot = tiledlayout(2,2,'TileSpacing','Compact');
set(gcf,'Position',[300 100 800 600])
for aa = 1:4
    fval = data{aa};
    xXaxis = fval.funval(:,1);
    yYaxis = fval.funval(:,2);
    ax(aa)=nexttile(aa);
    set(ax(aa),...
        'box','on',...
        'FontSize',fontsize);
    set(gca,'TickLength',[0.025, 0.01])
    hold on

    plot(ax(aa),xXaxis,yYaxis,'Color',[0 0 0 1],'LineWidth',1,'LineStyle', '-')
    plot(ax(aa),xXaxis,threshold*ones(size(xXaxis)),'Color',[0 0 0 1],'LineStyle','--','LineWidth',1)
    plot(ax(aa),MLpars(aa),-MLvalue,...
        'Marker','o',...
        'MarkerSize',8,...
        'MarkerFaceColor',[1 0 0],...
        'MarkerEdgeColor',[1 0 0],...
        'LineStyle','none')
    hold on

    xlabel(ax(aa),xLabels{aa},'FontSize',14,'Interpreter','latex');
    ylabel(ax(aa),'$-$ Profile likelihood','FontSize',14,'Interpreter','latex');
    text(0.075,0.9,figlabs{aa},...
        'Units','Normalized',...
        'HorizontalAlignment','left',...
        'FontSize',fontsize,...
        'FontWeight','Normal'); 

end

tiledplot.TileSpacing = 'compact';
tiledplot.Padding = 'compact';

savefig('./figures/practical_identifiability.fig')
exportgraphics(gcf,...
    './figures/practical_identifiability.png')
exportgraphics(gcf,...
    './figures/practical_identifiability.eps',...
    'ContentType','vector')



end