data=get_ifndata(14)
data_54=data(1:7,2:3);
data_31=data(8:11,2:3);
data_14=data(12:17,2:3);
data_11=data(18:23,2:3);
data_8=data(24:28,2:3);
%plot(data_54(:,1),data_54(:,2))
%hold on
%plot(data_31(:,1),data_31(:,2))
%plot(data_14(:,1),data_14(:,2))
%plot(data_11(:,1),data_11(:,2))
%plot(data_8(:,1),data_8(:,2))
%hold off
%Hill(x,y,noParam,maximum,slope,halfActiv,intercept);
%Hill(data_54(:,1), data_54(:,2),3,3,10,2.5,0)
x_54=(data_54(:,1))
x_31=(data_31(:,1))
x_14=(data_14(:,1))
x_11=(data_11(:,1))
x_8=(data_8(:,1))
% MAPPING: Emax = b(1),  EC50 = b(2)
err=[];
for qc = 1:10
    hill_fit = @(b,x)  b(1).*(x.^qc)./(b(2)^2+(x.^qc));
    b0 = [4; 9];                                  % Initial Parameter Estimates
    B = lsqcurvefit(hill_fit, b0, x_54 ,data_54(:,2));
    AgVct = linspace(min(x_54), max(x_54 ));   % Plot Finer Resolution
    B_31 = lsqcurvefit(hill_fit, b0, x_31 ,data_31(:,2));
    AgVct_31 = linspace(min(x_31 ), max(x_31 ));   % Plot Finer Resolution
    B_14 = lsqcurvefit(hill_fit, b0, x_14 ,data_14(:,2));
    AgVct_14 = linspace(min(x_14 ), max(x_14 ));
    B_11 = lsqcurvefit(hill_fit, b0, x_11 ,data_11(:,2));
    AgVct_11 = linspace(min(x_11 ), max(x_11 ));
    B_8 = lsqcurvefit(hill_fit, b0,  x_8 ,data_8(:,2));
    AgVct_8 = linspace(min(x_8) , max(x_8));
    di_54=data_54(:,2)-hill_fit(B,x_54);
     di_31=data_31(:,2)-hill_fit(B_31,x_31);
      di_14=data_14(:,2)-hill_fit(B_14,x_14);
       di_11=data_11(:,2)-hill_fit(B_11,x_11);
        di_8=data_8(:,2)-hill_fit(B_8,x_8);

    err(qc)=norm(di_54)+norm(di_31)+norm(di_14)+norm(di_11)+norm(di_8)

end
hill_fit = @(b,x)  b(1).*(x.^2)./(b(2)^2+(x.^2));
    b0 = [4; 9];                                  % Initial Parameter Estimates
    B = lsqcurvefit(hill_fit, b0, x_54 ,data_54(:,2));
    AgVct = linspace(min(x_54), max(x_54 ));   % Plot Finer Resolution
    B_31 = lsqcurvefit(hill_fit, b0, x_31 ,data_31(:,2));
    AgVct_31 = linspace(min(x_31 ), max(x_31 ));   % Plot Finer Resolution
    B_14 = lsqcurvefit(hill_fit, b0, x_14 ,data_14(:,2));
    AgVct_14 = linspace(min(x_14 ), max(x_14 ));
    B_11 = lsqcurvefit(hill_fit, b0, x_11 ,data_11(:,2));
    AgVct_11 = linspace(min(x_11 ), max(x_11 ));
    B_8 = lsqcurvefit(hill_fit, b0,  x_8 ,data_8(:,2));
    AgVct_8 = linspace(min(x_8) , max(x_8));

figure()
plot(AgVct, hill_fit(B,AgVct),'-r')
hold on

legend('54')
plot(x_54 ,data_54(:,2), 'rp','HandleVisibility','off')
plot(x_31 ,data_31(:,2), 'bp','HandleVisibility','off')
plot(x_14(:,1),data_14(:,2), 'yp','HandleVisibility','off')
plot(x_11(:,1),data_11(:,2), 'gp','HandleVisibility','off')
plot(x_8(:,1),data_8(:,2), 'mp','HandleVisibility','off')
plot(AgVct_31, hill_fit(B_31,AgVct_31), '-b','DisplayName','31')

plot(AgVct_14, hill_fit(B_14,AgVct_14), '-y','DisplayName','14')

plot(AgVct_11, hill_fit(B_11,AgVct_11), '-g','DisplayName','11')

plot(AgVct_8, hill_fit(B_8,AgVct_8), '-m','DisplayName','8' )
 
grid
xlabel('nt')
ylabel('infg')

hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
hill__11 = @(x)  B_11(1).*(x.^2)./((B_11(2))^2+(x.^2));
hill__8 = @(x)  B_8(1).*(x.^2)./((B_8(2))^2+(x.^2));
hill__54= @(x)  B(1).*(x.^2)./((B(2))^2+(x.^2));
% -------------------------------------------------------------------------
% Algorithm to fit data by estimating parameters with MCMC.
% -------------------------------------------------------------------------
%
% This code estimates the parameters for the model through an MCMC
% simulation. There are 3 possible sets of parameters to fit:
% (1) KD, kp, kappa;
% (2) KD, kp, kappa, mu;
% (3) KD, pi, mu, alpha, C;
% The model fitted here is averaged (i.e. does not consider the
% distributions).
% The MCMC is initially run with the set of parameters derived from
% estimations of the activation function. It is afterwards re-run with the
% parameters that generated the best fit in the first MCMC estimation.
%
% In the model considered here, we account for the number of engaged TCR,
% keeping record of serial engagement.
% 
% The data corresponding to both radii is simultaneously fit, to make sure
% that all parameters are consistent throughout the experiments.
% -------------------------------------------------------------------------

clear data model params options

% Cluster Paths
addpath(['../project_code/mcmcstat-master']);
addpath('ifngProd','dist_fitting','dyn_models','dyn_models_NP','zeroNewton','insProb_functions');

%data = [[repmat(14,28,1), get_ifndata(14)];[repmat(20,21,1), get_ifndata(20)]];
%data20 = [[repmat(20,21,1), get_ifndata(20)]];
%data14 = [repmat(14,28,1), get_ifndata(14)];
%data14(1:17,:) = [];
%data = data14;
data= rand(300,4);
data(:,1)=14;

data(1:100,2)=8;
data(101:200,2)=11;
data(201:300,2)=54;
data(1:100,3)=linspace(10^-1,50,100)
data(101:200,3)=linspace(10^-1,50,100)
data(201:300,3)=linspace(10^-1,50,100)
data(:,4)=[hill__8(data(1:100,3));hill__11(data(101:200,3));hill__54(data(201:300,3))]




%% Initial Parameter Values

%param0 = [19, 0.5, 9.9, 0.18, 4.7]; %[19, 0.5, 9.9, 0.18, 4.7, 5.0]; % [kd, {pt; KP, N}, mu, alpha, C, (angle)]
param0 = [0.3487,    8.9550  ,  0.1050   ,47.5631  , 20.6581];

%param0 = [19, 0.5, 9.9, 0.18, 4.7, 1]; % (kd, pt, mu, alpha, C, KD) where KD is CA passage ratio.

% Computation of initial variance
ss0 = louis_errorFunction(param0,data);

%% Settings for MCMC
params = {
    {'kd',param0(1),eps,1}
    
    %{'pt',param0(2),eps,1}
    %{'phi',param0(1),eps,3000}
    %{'b',param0(2),eps,3000}
    
    {'mu',param0(2),5,20}                % (5, 20)
    {'alph',param0(3),eps,1}            % (eps,10); If using hill function set alpha = (eps,5)
    {'C',param0(4),eps,1e3}               % (eps,20)
    {'KD',param0(5),eps,1e3}
    };

% Error settings
model.sigma2 = ss0;%/(length(data.ydata(:,1))-1);
model.ssfun = @louis_errorFunction;

% Simulation settings
nsim=150; %1500;
options.nsimu = nsim;
options.updatesigma = 1;
options.method = 'dram';
%options.savesize = 500;
%options.savedir = './TEMP_MCMC/';
%% Run MCMC algorithm
[result, chain, s2chain, sschain] = mcmcrun(model,data,params,options);

%% Save MCMC algorithm results
mkdir([pp,'FUNCfig_nsim',num2str(nsim),'date',date,'/'])
filename = [pp,'FUNCfig_nsim',num2str(nsim),'date',date,'/results.mat'];
%save(filename,'result','chain','s2chain','sschain')
%--------------------------------------------------------------------------
% For Saving on the Cluster (must save in --ascii format)

chain_name = [pp,'fig_nsim',num2str(nsim),'date',date,'/chain','.txt'];
ss_name = [pp,'fig_nsim',num2str(nsim),'date',date,'/sschain','.txt'];
s2_name = [pp,'fig_nsim',num2str(nsim),'date',date,'/s2chain','.txt'];

save(chain_name,'chain','-ascii')
save(ss_name,'sschain','-ascii')
save(s2_name,'sschain','-ascii')
