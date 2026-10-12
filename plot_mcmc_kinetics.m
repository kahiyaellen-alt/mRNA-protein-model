function [] = plot_mcmc_kinetics()

clc
close all
set(0,'DefaultFigureVisible','on');

%% Load synthetic data
metadata = load("./synthetic_data.mat");

time_data = metadata.time_data;
time_curve = metadata.time_curve;
mRNA_synthetic = metadata.mRNA_synthetic;
protein_synthetic = metadata.protein_synthetic;
tmax = 24;

%% Initial conditions
y0 = [2.5; 6.5];

% %% Load MCMC chain
 chainsStruct = load('./chains/chains.mat');
 chains = chainsStruct.chains;
 
 %% Thin chain
 keep = 1:2:size(chains,1);
 nSamples = length(keep);

%% Store posterior trajectories
 mRNA_all    = zeros(nSamples,length(time_curve));
 protein_all = zeros(nSamples,length(time_curve));

 for ii = 1:nSamples
 
     pars = chains(keep(ii),:);
 
     output = model(pars);
 
     mRNA_all(ii,:)    = output{1};
     protein_all(ii,:) = output{2};
 
end

%% Min-Max envelopes
mRNA_min = min(mRNA_all,[],1);
mRNA_max = max(mRNA_all,[],1);

protein_min = min(protein_all,[],1);
protein_max = max(protein_all,[],1);

%% True parameter values
true_pars = [20 0.5 0.25 0.5];

output_true = model(true_pars);

mRNA_true    = output_true{1};
protein_true = output_true{2};

%% Plot
%% ============================================================

%figure
tiledplot = tiledlayout(1,2,'TileSpacing','Compact');
set(gcf,'Position',[200 200 800 300])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mRNA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax(1) = nexttile(1);
hold(ax(1),'on')
box(ax(1),'on')
patch(ax(1),...
    [time_curve fliplr(time_curve)],...
    [mRNA_min fliplr(mRNA_max)],...
    [0 0.4470 0.7410],...
    'FaceAlpha',0.2,...
    'EdgeColor','none');
plot(ax(1),...
    time_curve,...
    mRNA_true,...
    'k-',...
    'LineWidth',2)
plot(ax(1),...
    time_data,...
    mRNA_synthetic,...
    'ro',...
    'MarkerSize',8,...
    'MarkerFaceColor','r',...
    'LineWidth',1.5)
xlabel(ax(1),'Time (h)')
ylabel(ax(1),'mRNA')
title(ax(1),'mRNA')
text(0.02, 0.95, 'A1', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Protein
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax(2) = nexttile(2);
hold(ax(2),'on')
box(ax(2),'on')
patch(ax(2),...
    [time_curve fliplr(time_curve)],...
    [protein_min fliplr(protein_max)],...
    [0 0.4470 0.7410],...
    'FaceAlpha',0.2,...
    'EdgeColor','none');
plot(ax(2),...
    time_curve,...
    protein_true,...
    'k-',...
    'LineWidth',2)
plot(ax(2),...
    time_data,...
    protein_synthetic,...
    'ro',...
    'MarkerSize',8,...
    'MarkerFaceColor','r',...
    'LineWidth',1.5)
xlabel(ax(2),'Time (h)')
ylabel(ax(2),'Protein')
title(ax(2),'Protein')
text(0.02, 0.95, 'A2', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shared title
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tiledplot.TileSpacing = 'compact';
tiledplot.Padding = 'compact';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('./figures','dir')
    mkdir('./figures')
end

savefig('./figures/kinetics.fig')
exportgraphics(gcf,'./figures/kinetics.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve ODE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function output = model(pars)

        sol = ode45(@(t,y) des(t,y,pars),...
                    [0 tmax],...
                    y0);

        output{1} = deval(sol,time_curve,1);
        output{2} = deval(sol,time_curve,2);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = des(~,y,pars)

        mRNA    = y(1);
        protein = y(2);

        dydt = [pars(1) - pars(2)*mRNA;...
                pars(3)*mRNA - pars(4)*protein];

    end

end