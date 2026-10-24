function [] = cornerplot()

%%% command:
%%% cornerplot()

clc;
close all;
set(0,'DefaultFigureVisible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD MCMC SAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load the saved MCMC chain
chainData = load('./chains/chains.mat');
chains = chainData.chains;

chains = log10(chains);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER LABELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a{1} = 'k_1';
a{2} = 'k_2';
a{3} = 'k_3';
a{4} = 'k_4';

labs = {a{1}, ...
    a{2}, ...
    a{3}, ...
    a{4}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATE CORNER PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Requires the gwmcmc package with ecornerplot
figure();

ecornerplot(chains, ...
    'ks',true, ...
    'color',[.3 .3 .3], ...
    'names',labs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('./figures','dir')
    mkdir('./figures');
end

savefig('figures/corner.fig');
exportgraphics(gcf,'figures/corner.png');

end
