function []= synthetic_data()

clc;
close all;

format long

%% Time points in hrs
tmax = 24;
time_data = [0:0.5:5,6:4:tmax];
time_curve = 0:0.1:tmax;

%% Initial conditions
y0 = [2.5; 6.5];

%% True parameters
k1 = 20;
k2 = 0.5;
k3 = 0.25;
k4 = 0.5;

true_pars = [k1 k2 k3 k4];

%% Steady states for fold-change
Mstar = 40;   % mRNA steady state
Pstar = 20;   % protein steady state

%% Solve ODE system
sol = ode45(@(t,y) model(t,y,true_pars), [0 tmax], y0);

%% Solve ODE system for time curve 
[T,Y] = ode45(@(t,y) model(t,y,true_pars), time_curve, y0);
mRNA_curve = Y(:,1);
protein_curve = Y(:,2);

%% True model outputs evaluated at measured time points
mRNA_data    = deval(sol,time_data,1);
protein_data = deval(sol,time_data,2);

%% Gaussian parameters
mu    = 0;
stdev_mRNA = 2.5;
stdev_protein = 1;

%% Generate synthetic noisy data
mRNA_noise    = normrnd(mu, stdev_mRNA, size(mRNA_data));
protein_noise = normrnd(mu, stdev_protein, size(protein_data));

%% Add noise to model output
mRNA_synthetic    = mRNA_data + mRNA_noise;
protein_synthetic = protein_data + protein_noise;

%% Ensure no negative values
mRNA_synthetic(mRNA_synthetic < 0) = 0;
protein_synthetic(protein_synthetic < 0) = 0;

%% Convert to fold-change relative to steady states
mRNA_fc_curve    = mRNA_curve / Mstar;
protein_fc_curve = protein_curve / Pstar;

mRNA_fc_data    = mRNA_synthetic / Mstar;
protein_fc_data = protein_synthetic / Pstar;

%% Display synthetic fold-change data
disp('Synthetic mRNA fold-change data:')
disp(mRNA_fc_data)

disp('Synthetic protein fold-change data:')
disp(protein_fc_data)

% Save synthetic data to .mat file 
outname = ['synthetic_data'  '.mat'];  

% Variables to save
saveVars = { ...
    'time_data', ...
    'mRNA_synthetic', ...
    'protein_synthetic', ...
    'mRNA_fc_data', ...
    'protein_fc_data', ...
    'time_curve', ...
    'mRNA_curve', ...
    'protein_curve', ...
    'mRNA_fc_curve', ...
    'protein_fc_curve', ...
    'Mstar', ...
    'Pstar' ...
     };

% Save to file
save(outname, saveVars{:});

fprintf('Saved synthetic data to %s\n', outname);


%% Plot results (fold-change)
figure

subplot(2,1,1)
plot(T, mRNA_fc_curve, 'b-', 'LineWidth',2)
hold on
plot(time_data, mRNA_fc_data, 'ro', 'MarkerSize',8,'LineWidth',1.5)
xlabel('Time')
ylabel('mRNA (fold-change)')
title('Synthetic mRNA Fold-Change Data')
legend('True model (fold-change)','Synthetic data (fold-change)')

subplot(2,1,2)
plot(T, protein_fc_curve, 'b-', 'LineWidth',2)
hold on
plot(time_data, protein_fc_data, 'ro', 'MarkerSize',8,'LineWidth',1.5)
xlabel('Time')
ylabel('Protein (fold-change)')
title('Synthetic Protein Fold-Change Data')
legend('True model (fold-change)','Synthetic data (fold-change)')

%%% save figure %%%
if ~exist('./figures', 'dir')
    mkdir('./figures')
end
savefig(strcat('figures/synthetic_data.fig'));
exportgraphics(gcf,'figures/synthetic_data.png');

%% ODE Model 
function dydt = model(~,y,pars)

mRNA   = y(1);
protein = y(2);

k1 = pars(1);
k2 = pars(2);
k3 = pars(3);
k4 = pars(4);

dydt = [k1 - k2*mRNA;...
    k3*mRNA - k4*protein];

end

end
