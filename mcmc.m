function [] = mcmc()

 metadata = load("./synthetic_data.mat");

    data    = {metadata.mRNA_synthetic;...
               metadata.protein_synthetic};
    stdev   = {2.5;...
               1};
    tmax    = 24;
    y0 = [2.5; 6.5];
    time    = metadata.time_data;

    k1 = 20;
    k2 = 0.5;
    k3 = 0.25;
    k4 = 0.5;

    pars = [k1 k2 k3 k4];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SLICE SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        nSamples = 50000;
        burnIn = 2000;
        thin = 5;
        chains = slicesample(pars,nSamples,'logpdf',@logprob);
        chains = chains(burnIn+1:thin:end,:);
        
        
        %%% save chains %%%
        if ~exist('./chains', 'dir')
            mkdir('./chains')
        end

        save(strcat('./chains/chains.mat'),'chains');
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% helper functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = -0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev);
    end

    %%% log-probability %%%
    function [value] = logprob(pars)  
        lp = logprior(pars);
        if ~isfinite(lp)
            value = -inf;
        else
            value = lp+loglike(pars);
        end   
    end

    %%% priors %%%
    function [flag] = logprior(pars)
        if all(pars>0)
            flag = 0;
        else
            flag = -inf;
        end    
    end

    %%% log-likelihood %%%
    function [value] = loglike(pars)  
        [solution] = model(pars);
        value = 0;
        for ii=1:length(data)
              value = value + sum(lognormpdf(data{ii}, solution{ii}, stdev{ii}));      
        end
    end

    %%% output of the model %%%
    function [output] = model(pars)
        solution = ode45(@(t,y) des(t,y,pars), [0 tmax], y0);

        output{1} = deval(solution,time,1);
        output{2} = deval(solution,time,2);
        
    end


    function dydt = des(~,y,pars)
        
        zeros(2,1);

        mRNA   = y(1);
        protein = y(2);

        

        dydt = [pars(1) - pars(2)*mRNA;...
                pars(3)*mRNA - pars(4)*protein];


    end
    
    
end 