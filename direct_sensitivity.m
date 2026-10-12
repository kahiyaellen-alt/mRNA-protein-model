function [] = direct_sensitivity()

%%% command:
%%% direct_sensitivity()  

clc;
close all
set(0,'DefaultFigureVisible','on');

% parameters
k_1 = 20;
k_2 = 0.5;
k_3 = 0.25;
k_4 = 0.5;

pars = [k_1 k_2 k_3 k_4];

parnames{1} = '$k_\mathrm{1}$';
parnames{2} = '$k_\mathrm{2}$';
parnames{3} = '$k_\mathrm{3}$';
parnames{4} = '$k_\mathrm{4}$';

tmax    = 24;

Text  = {'$mRNA$';...
    '$Protein$'};

if ~exist('./figures', 'dir')
    mkdir('./figures')
end

%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
colors = [0 0 0; 1 0 0; 0 0 1; 0 1 0];

% Use parnames directly for legend entries (they are strings already)
for ii = 1:length(pars)
    Legend{ii} = parnames{ii};
end

figlabs = {'(A)','(B)'};
tiledplot = tiledlayout(1,2,'TileSpacing','Compact');
set(gcf,'Position',[100 100 800 300])

mm = 0;
for aa = 1:2
    mm = mm+1;
    ax(mm) = nexttile(mm);
    set(ax(mm),...
        'box','on',...
        'FontSize',10);
    set(gca,'TickLength',[0.025, 0.01])

    hold on;
    for i = 1:length(parnames)
        [tt, output] = model(pars,i);
        plot(ax(mm),tt,output(:,mm),'Color',colors(i,:),'LineWidth',1.5,'LineStyle','-')
        hold on;
    end
    title(strcat('\rm',Text{mm}),'FontSize',10,'Interpreter','latex');
    text(0.95, 0.96,figlabs{mm},...
        'Units','Normalized',...
        'HorizontalAlignment','right',...
        'FontSize',10,...
        'FontWeight','Normal');
    xlabel(ax(mm),'Time (hours)','FontSize',10);
end

leg = legend(ax(mm),Legend{:,:},...
    'Location','Eastoutside',...
    'FontSize',10,...
    'Orientation','Vertical',...
    'Interpreter','latex');
% rect = [0.5, 0.6, .15, .2];
% set(leg, 'Position', rect)
set(leg, 'Box', 'off');

title(leg,'Parameter')
ylabel(tiledplot,'Normalized sensitivity coefficient','FontSize',10);
tiledplot.TileSpacing = 'compact';
tiledplot.Padding = 'compact';

savefig('figures/direct_sensitivity.fig');
exportgraphics(gcf,'figures/direct_sensitivity.png');
exportgraphics(gcf,'figures/direct_sensitivity.eps','ContentType','vector');


    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% output of the model %%%
    function [t,output] = model(pars,parsnum) 
        [t,s] = ode23s(@ode,[0 tmax],...
                       [2.5 6.5 0 0],...
                       [],...
                       pars,...
                       parsnum); 
        output_m = s(:,1:2);
        output_s = s(:,3:end);
        output = pars(parsnum) * output_s ./ output_m;
    end

    %%% Jacobian times sensitivities %%%
    function [dsdt] = ode(~,x,pars,parsnum) 

        dfdp_matrix = [1 -x(1) 0 0;...
                        0 0 x(1) -x(2)];
        dfdp = [0;0;dfdp_matrix(:,parsnum)];
        
        dsdt = [pars(1)- pars(2)*x(1);...                     
                pars(3)*x(1)- pars(4)*x(2);...                                                            
                -pars(2)*x(3);...
                pars(3)*x(3)-pars(4)*x(4)] + dfdp;       
    end

end