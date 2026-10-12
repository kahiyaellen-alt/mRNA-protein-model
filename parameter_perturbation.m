function [] = parameter_perturbation()

    %%% Command:
    %%% parameter_perturbation()  
   
    clc;
    close all;
    set(0, 'DefaultFigureVisible', 'on');

    metadata = load("./synthetic_data.mat");

    % parameters
    k_1 = 20;
    k_2 = 0.5;
    k_3 = 0.25;
    k_4 = 0.5;

    pars = [k_1 k_2 k_3 k_4];
    tmax    = 24;
    T = metadata.time_curve;

    Text = {'mRNA'; 'Protein'; 
            'mRNA'; 'Protein';
            'mRNA'; 'Protein';
            'mRNA'; 'Protein'};

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
param_names = {'$k_\mathrm{1}$',...
                '$k_\mathrm{1}$',...
                '$k_\mathrm{2}$',...
                '$k_\mathrm{2}$',...
                '$k_\mathrm{3}$',...
                '$k_\mathrm{3}$',...
                '$k_\mathrm{4}$',...
                '$k_\mathrm{4}$'};

    %%% Sensitivity Analysis %%%
    sensitivity_range = [-0.5, 0, 0.5]; % -50%, best-fit, +50% 
    colors = [0 0 1; 0 0 0; 0 1 1]; % lower bound, best-fit, upper bound
    
    %%% Create one figure for all parameters
    tiledplot = tiledlayout(4,2,'TileSpacing','Compact');
    figlabs = {'(A)','(B)',...
               '(C)','(D)',...
               '(E)','(F)',...
               '(G)','(H)'};
    set(gcf, 'Position', [300, 100, 1000, 850]);
    num_rows = length(pars);
    num_cols = 2; 
    
    mm = 0;
    for param_idx = 1:num_rows
        for var_idx = 1:num_cols
            mm = mm+1;
            ax(mm) = nexttile(mm);
            set(ax(mm),...
                'box','on',...
                'FontSize',12);
            set(gca,'TickLength',[0.025, 0.01]);
            hold on;
        
            % Results storage
            outputs = cell(size(sensitivity_range));

            % Simulate for -50%, best-fit, and +50%  
            for i = 1:length(sensitivity_range)
                modified_pars = pars;
                modified_pars(param_idx) = pars(param_idx) * (1 + sensitivity_range(i));
                output = model(modified_pars, T, tmax);
                outputs{i} = output{var_idx};
            end

            % Plot curves
            lower_bound = outputs{1}; % -50%
            best_fit   = outputs{2}; % Best-fit
            upper_bound = outputs{3}; % +50%

            % Plot lower bound
            plot(T, lower_bound, '-','Color', colors(1, :), 'LineWidth', 1, 'DisplayName', '-50%');
            hold on;
            % Plot best-fit
            plot(T, best_fit, 'Color', colors(2, :), 'LineWidth', 1, 'DisplayName', 'best_fit');
            hold on;
            % Plot upper bound
            plot(T, upper_bound, '-', 'Color', colors(3, :), 'LineWidth', 1, 'DisplayName', '+50%');
            hold on;

            xlabel(ax(mm), 'Time (hours)');
           %title(ax(mm), strcat(Text{mm}), 'FontSize', 12, 'FontWeight', 'Normal');
           title(ax(mm), sprintf('%s, %s = %.3g', Text{mm}, param_names{mm}, pars(param_idx)), ...
               'FontSize', 12, 'FontWeight', 'Normal', 'Interpreter', 'latex');

            hold on;

            text(0.1,0.9,figlabs{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',12,...
                'FontWeight','Normal');    
            hold on;

            %title(ax(mm), Text{mm}, 'FontSize', 12, 'FontWeight', 'Normal');

            text(ax(mm), 0.95, 0.90, param_names{mm}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'right', ...
                'FontSize', 12, ...
                'FontWeight', 'Normal', ...
                'Interpreter', 'latex');
            

            hold on;
            legend('-50%', 'best-fit', '+50%', 'Location', 'northeast', 'FontSize', 10);
            % end
        end
    end
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    % Save the combined figure
    savefig('figures/combined_parameter_perturbation.fig');
    exportgraphics(gcf, 'figures/combined_parameter_perturbation.png');
    %exportgraphics(gcf, 'figures/combined_parameter_perturbation.eps', 'ContentType', 'vector');


    %%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [output] = model(pars, T, tmax)
        solution = ode23s(@ode, [0 tmax], [2.5 6.5], [], pars);
        for kk = 1:2
            output{kk} = deval(solution, T, kk);
        end
    end
    
    % ODE system
    function [dxdt] = ode(~, x, pars)
        dxdt = zeros(2, 1);

        dxdt(1) = pars(1)- pars(2)*x(1); % mRNA
        dxdt(2) = pars(3)*x(1)- pars(4)*x(2); % protein
    end

end