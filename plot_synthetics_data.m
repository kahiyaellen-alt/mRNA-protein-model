function []= plot_synthetics_data()

clc;
close all;

format long
metadata = load("./synthetic_data.mat");

%% Plot results
%figure
tiledplot = tiledlayout(1,2,'TileSpacing','Compact');
set(gcf,'Position',[200 200 800 300])
ax(1) = nexttile(1);
hold(ax(1),'on')
box(ax(1),'on')
plot(ax(1),metadata.time_curve, metadata.mRNA_curve, 'b-', 'LineWidth',2)
hold on
plot(ax(1),metadata.time_data, metadata.mRNA_synthetic, 'ro', 'MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.5)
xlabel('Time (h)')
ylabel('mRNA')
title('Synthetic mRNA Data')
legend('True model','Synthetic data')
text(0.02, 0.95, 'A1', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12)

ax(2) = nexttile(2);
hold(ax(2),'on')
box(ax(2),'on')
plot(metadata.time_curve,metadata.protein_curve, 'b-', 'LineWidth',2)
hold on
plot(metadata.time_data, metadata.protein_synthetic, 'ro', 'MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.5)
xlabel('Time (h)')
ylabel('Protein')
title('Synthetic Protein Data')
legend('True model','Synthetic data')
text(0.02, 0.95, 'A2', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12)

tiledplot.TileSpacing = 'compact';
tiledplot.Padding = 'compact';

%%% save figure %%%
if ~exist('./figures', 'dir')
    mkdir('./figures')
end
savefig(strcat('figures/synthetics_data.fig'));
exportgraphics(gcf,'figures/synthetics_data.png');
% exportgraphics(gcf,'../../LaTeX/figures/data.eps','ContentType','vector');



end