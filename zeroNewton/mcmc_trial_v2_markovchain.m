% -------------------------------------------------------------------------
% Algorithm to fit data by estimating parameters with MCMC.
% (this one uses the second version of binding)
% -------------------------------------------------------------------------
%% r = 14
clear data model params options
r = 20;
data = [[repmat(14,28,1), get_ifndata(14)];[repmat(20,21,1), get_ifndata(20)]];

% Parameters to fit (3 options):
% (1) KD, kp, kappa;
% (2) KD, kp, kappa, mu;
% (3) KD, pi, mu, alpha, C;
% % (1)
% param0 = [155,10,3];
% % (2)
% param0 = [55,100,3,8.08];
% % (3)
% param0 = [55,.95,8.08,0.39,5.93];
% param0 = [5,.24,14,0.39,5.93];
param0 = [23.131      0.51625       10.746      0.25838       5.2968];
ss0 = errorFunctionParFitting_v2_markovchain(param0,data)

% (1)
% params = {
%     {'KD',param0(1),eps,3000}
%     {'kp',param0(2),eps,10000}
%     {'kappa',param0(3),1,20}
%     };
% % (2)
% params = {
%     {'KD',param0(1),eps,3000}
%     {'kp',param0(2),eps,10000}
%     {'kappa',param0(3),1,20}
%     {'mu',param0(4),7.29,8.83}
%     };
% % (3)
params = {
    {'KD',param0(1),eps,3000}
    {'pi',param0(2),eps,1}
    {'mu',param0(3),7.29,15}%8.83}
    {'alpha',param0(4),0.17,0.6}
    {'C',param0(5),5.14,6.72}
    };

model.sigma2 = ss0;%/(length(data.ydata(:,1))-1);
model.ssfun = @errorFunctionParFitting_v2_markovchain;

nsim=3000;
options.nsimu = nsim;
options.updatesigma = 1;

[result, chain, s2chain, sschain] = mcmcrun(model,data,params,options);
%%
% filename = ['mcmcres',num2str(r),'_n',num2str(nsim),'multi.mat'];
% save(filename,'result','chain','s2chain','sschain');

f1 = figure(101);
mcmcplot(chain(round(nsim/2):end,:),[],result.names,'denspanel');
% saveas(f1,['denspanel',num2str(r),'_',num2str(nsim)],'jpg')

f2 = figure(102);
mcmcplot(chain(round(nsim/2):end,:),[],result.names,'pairs');
% saveas(f2,['pairs',num2str(r),'_',num2str(nsim)],'jpg')

f3 = figure(103);
mcmcplot(chain(round(nsim/2):end,:),[],result.names,'hist');
% saveas(f3,['hist',num2str(r),'_',num2str(nsim)],'jpg')

f4 = figure(104);
mcmcplot(chain(round(nsim/2):end,:),[],result.names,'chainpanel');
% saveas(f3,['hist',num2str(r),'_',num2str(nsim)],'jpg')
%%
minind = round(nsim/2); % startpoint of convergence
chainconv = chain(minind:end,:);
lconv = length(chainconv);

%randomized sampling of values
nsample = 20;
nindrnd = unifrnd(0,ones(nsample,1));
nind = floor(lconv*nindrnd);
chainsample = chainconv(nind,:);
paramsample = mean(chainsample);
predinfsample = predinf_vals(paramsample,data);
figure(1001)
plot_predinf(data,predinfsample);

%plot predictions for minimum error found
indmin = find(sschain==min(sschain));
parammin = chain(indmin);
predinfmin = predinf_vals(parammin,data);
figure(1002)
plot_predinf()

%nsample = 500;
% 
% 
% modelfun = @predinf;
% out = mcmcpred(result,chain,s2chain,data.ydata,modelfun,nsample);
% 



%% r = 20

r = 20;
data2.ydata = get_ifndata(r);

predinf2 = predinf_vals([1e-3,1e-2,0.55,10000],data2)
[param02, ss02] = fminsearch(@(param)predinf_err(param,data),[1e-3,1e-2,.1,3000])

