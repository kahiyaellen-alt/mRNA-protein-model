function [] = practical_identifiability()

clc;
clear all;
close all;


%% Load data
metadata = load("./synthetic_data.mat");

data = {metadata.mRNA_synthetic;...
        metadata.protein_synthetic};

stdev = {2.5;...
         1};

tmax = 24;

y0 = [2.5;6.5];

time = metadata.time_data;


%% Maximum likelihood parameter values from MCMC
chains = load('./chains/chains.mat');
chains = chains.chains;

MLvalue = -10000;
MLpars = [1, 1, 1, 1];
dim = size(chains,1);
for i =1:dim
    ipars = chains(i,:);
    loglikevalue = loglike(ipars);
    if loglikevalue > MLvalue
        MLvalue = loglikevalue;
        MLpars = ipars;
    end
end
disp(MLvalue);
disp(MLpars);

save('./MLvalue.mat','MLvalue')
save('./MLpars.mat','MLpars')
 
% %% Parameter ranges for profile likelihood
k_1_span = 15:0.5:25;
k_2_span = 0.2:0.02:1.0;
k_3_span = 0.1:0.01:0.5;
k_4_span = 0.2:0.02:1.0;
 
%%% k_1 %%%
funval = [];    
options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point');
for i=1:length(k_1_span)
    pars = [MLpars(2), MLpars(3), MLpars(4)];
    [x,fval,exitflag,output] = fmincon(@(params)loglikek_1(params,k_1_span(i)), pars, [], [], [], [], 0*ones(length(pars),1), inf*ones(length(pars),1), [], options);
    funval = [funval; k_1_span(i), fval];
    disp(['k_1 ',num2str(i),' done. Value = ',num2str(fval)])
end    
save('./k_1.mat','funval');
disp('k_1 done.')

%%% k_2 %%%
funval = [];
for i=1:length(k_2_span)
    % params are [k1, k3, k4] (k2 fixed)
    pars = [MLpars(1), MLpars(3), MLpars(4)];
    [x,fval,exitflag,output] = fmincon(@(params)loglikek_2(params,k_2_span(i)), pars, [], [], [], [], 0*ones(length(pars),1), inf*ones(length(pars),1), [], options);
    funval = [funval; k_2_span(i), fval];
    disp(['k_2 ',num2str(i),' done. Value = ',num2str(fval)])
end
save('./k_2.mat','funval');
disp('k_2 done.')

%%% k_3 %%%
funval = [];
for i=1:length(k_3_span)
    % params are [k1, k2, k4] (k3 fixed)
    pars = [MLpars(1), MLpars(2), MLpars(4)];
    [x,fval,exitflag,output] = fmincon(@(params)loglikek_3(params,k_3_span(i)), pars, [], [], [], [], 0*ones(length(pars),1), inf*ones(length(pars),1), [], options);
    funval = [funval; k_3_span(i), fval];
    disp(['k_3 ',num2str(i),' done. Value = ',num2str(fval)])
end
save('./k_3.mat','funval');
disp('k_3 done.')

%%% k_4 %%%
funval = [];
for i=1:length(k_4_span)
    % params are [k1, k2, k3] (k4 fixed)
    pars = [MLpars(1), MLpars(2), MLpars(3)];
    [x,fval,exitflag,output] = fmincon(@(params)loglikek_4(params,k_4_span(i)), pars, [], [], [], [], 0*ones(length(pars),1), inf*ones(length(pars),1), [], options);
    funval = [funval; k_4_span(i), fval];
    disp(['k_4 ',num2str(i),' done. Value = ',num2str(fval)])
end
save('./k_4.mat','funval');
disp('k_4 done.')

%%log-likelihood for k_1
function [value] = loglikek_1(pars,k_1)
    pars_k1=[k_1,pars];
    solution=model(pars_k1);
    value=0;
    for ii=1:length(data)
        value=value+sum(lognormpdf(data{ii},solution{ii},stdev{ii}));
    end
       value=-value;
end

%%log-likelihood for k_2
function [value] = loglikek_2(pars,k_2)
    % pars = [k1, k3, k4]
    pars_k2 = [pars(1), k_2, pars(2), pars(3)];
    solution = model(pars_k2);
    value = 0;
    for ii=1:length(data)
        value = value + sum(lognormpdf(data{ii}, solution{ii}, stdev{ii}));
    end
    value = -value;
end

%%log-likelihood for k_3
function [value] = loglikek_3(pars,k_3)
    % pars = [k1, k2, k4]
    pars_k3 = [pars(1), pars(2), k_3, pars(3)];
    solution = model(pars_k3);
    value = 0;
    for ii=1:length(data)
        value = value + sum(lognormpdf(data{ii}, solution{ii}, stdev{ii}));
    end
    value = -value;
end

%%log-likelihood for k_4
function [value] = loglikek_4(pars,k_4)
    % pars = [k1, k2, k3]
    pars_k4=[pars,k_4];
    solution = model(pars_k4);
    value = 0;
    for ii=1:length(data)
        value = value + sum(lognormpdf(data{ii}, solution{ii}, stdev{ii}));
    end
    value = -value;
end

%%% log-likelihood formula %%%
function [value] = lognormpdf(data,model,stdev)
 value = -0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev);
end
    
%%% log-likelihood %%%
function [value] = loglike(pars)  
    [solution] = model(pars);
    value = 0;
    for ii=1:length(data)
        value = value + sum(lognormpdf(data{ii}, solution{ii}, stdev{ii}));
    end
end   

% %%% Solve ODE model
function output=model(pars)
solution=ode45(@(t,y)ode(t,y,pars),...
               [0 tmax],...
               y0);
output{1}=deval(solution,time,1);
output{2}=deval(solution,time,2);
end

%%% ODE system
function dydt=ode(~,y,pars)

mRNA=y(1);
protein=y(2);

k1=pars(1);
k2=pars(2);
k3=pars(3);
k4=pars(4);

dydt=[k1-k2*mRNA;...
      k3*mRNA-k4*protein];


end






end
