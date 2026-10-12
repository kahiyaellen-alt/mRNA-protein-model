function []= plot_synthetic_data()

clc;
close all;

format long
metadata = load("./synthetic_data.mat");


%% Plot results
figure

subplot(2,1,1)
plot(metadata.time_curve, metadata.mRNA_curve, 'b-', 'LineWidth',2)
hold on
plot(metadata.time_data, metadata.mRNA_fc_data, 'ro', 'MarkerSize',8,'LineWidth',1.5)
xlabel('Time')
ylabel('mRNA')
title('Synthetic mRNA Data')
legend('True model','Synthetic data')

subplot(2,1,2)
plot(metadata.time_curve,metadata.protein_curve, 'b-', 'LineWidth',2)
hold on
plot(metadata.time_data, metadata.protein_fc_data, 'ro', 'MarkerSize',8,'LineWidth',1.5)
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



end