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

%% Solve ODE system
sol = ode45(@(t,y) model(t,y,true_pars), [0 tmax], y0);
%t = sol.x;
%z = sol.y;

%% Solve ODE system for time curve 
[T,Y] = ode45(@(t,y) model(t,y,true_pars), time_curve, y0);
mRNA_curve = Y(:,1);
protein_curve = Y(:,2);

%% True model outputs
%mRNA_true    = z(1,:);
%protein_true = z(2,:);

%% True model outputs evaluated at measured time points
mRNA_data    = deval(sol,time_data,1);
protein_data = deval(sol,time_data,2);

%% Gaussian parameters
mu    = 0;
stdev_mRNA = 2.5;
stdev_protein = 1;

%% Generate synthetic noisy data
% Noise follows N(0, 0.1^2)

mRNA_noise    = normrnd(mu, stdev_mRNA, size(mRNA_data));
protein_noise = normrnd(mu, stdev_protein, size(protein_data));

%% Add noise to model output
mRNA_synthetic    = mRNA_data + mRNA_noise;
protein_synthetic = protein_data + protein_noise;

%% Ensure no negative values
mRNA_synthetic(mRNA_synthetic < 0) = 0;
protein_synthetic(protein_synthetic < 0) = 0;

%% Display synthetic data
disp('Synthetic mRNA data:')
disp(mRNA_synthetic)

disp('Synthetic protein data:')
disp(protein_synthetic)

% Save synthetic data to .mat file 

outname = ['synthetic_data'  '.mat'];  

% Variables to save (pick what you need)
saveVars = { ...
    'time_data', ...
    'mRNA_synthetic', ...
    'protein_synthetic', ...
    'time_curve', ...
    'mRNA_curve', ...
    'protein_curve'...
     };

% Save to file
save(outname, saveVars{:});

fprintf('Saved synthetic data to %s\n', outname);


%% Plot results
figure

subplot(2,1,1)
plot(T, mRNA_curve, 'b-', 'LineWidth',2)
hold on
plot(time_data, mRNA_synthetic, 'ro', 'MarkerSize',8,'LineWidth',1.5)
xlabel('Time')
ylabel('mRNA')
title('Synthetic mRNA Data')
legend('True model','Synthetic data')

subplot(2,1,2)
plot(T, protein_curve, 'b-', 'LineWidth',2)
hold on
plot(time_data, protein_synthetic, 'ro', 'MarkerSize',8,'LineWidth',1.5)
xlabel('Time')
ylabel('Protein')
title('Synthetic Protein Data')
legend('True model','Synthetic data')

%%% save figure %%%
if ~exist('./figures', 'dir')
    mkdir('./figures')
end
savefig(strcat('figures/synthetic_data.fig'));
exportgraphics(gcf,'figures/synthetic_data.png');
% exportgraphics(gcf,'../../LaTeX/figures/data.eps','ContentType','vector');

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