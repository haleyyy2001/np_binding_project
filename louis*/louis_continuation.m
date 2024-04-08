


function [kdGrid,stateKdContGrid] = louis_continuation(r,v,D,theta_deg,pp)
%% Initialize Parameters
pp = ['louis_results/'];          %Path to insProb results
%pp = ['louis_temp\theta',num2str(theta_deg),'\'];

load([pp,'insprob_r',num2str(r),'_v',num2str(v),'.mat']);           %Insertion probability of TCR onto NP
load([pp,'insprobNP_r',num2str(r),'.mat']);                         %Insertion probaability of NP onto NC

M = 1;%2000                         %Total Number of Nano-clusters
Mvec = M * probBindNPi(:);          %Proportion of NC reaching certain state


%% Contact Area Parameters
thefta = theta_deg * pi / 180;       %Convert Angle to Radians
rc = r * sin(theta_deg);                %Radius of contact area on Tcell surface
SAeff = 2*pi*(r^2)*(1-cos(theta_deg));  %Effective surface area of NP 
SAtot = 4*pi*(r^2);                 %Total surface area of NP
Acont = pi * r ^ 2;                 %Contact area on Tcell Surface
vh = round(SAeff * v / SAtot);      %Effective valency of NP

%% Continuation Parameters

konInitial = 0.1; 
koffInitial = 0.1;                              %Initial values of kon, koff and kd
kdInitial = koffInitial / konInitial;

kdMin = 1e-2;
kdMax = 1e4;                                    %Newton algorithm parameters
contStep = 1.15;
newtonTol = 1e-6;

%% Find Steady State of Initial (kon,koff)

mcov = 1;
for i=1:length(distCov)
    mcov = max(mcov,length(distCov{i}));
end
temp = dynsys_means(konInitial,koffInitial,v,vh,40,insProbTCRtoNP);             %dynsis_means(kf,kb,v,vh,Tm,pia)
B = unbindingRate(distCov,temp(:,3));

R = konInitial * v * D;

[t,y] = ode15s(@(t,x)f(R,Mvec,B,x,t),[0,500],[M,zeros(1,length(insProbNPtoNC))]');%,opts); f.m is the ODE
initState = y(end,:);

% Run continuation with respect to parameter kon (1), koff (2) or kd (3)
% indicated in the first argument

[kdGrid, stateKdContGrid] = continuation(3,kdMin,kdMax,contStep,konInitial,koffInitial,...
    initState,newtonTol,D,v,vh,insProbTCRtoNP,distCov,Mvec);

if exist([pp,'Continuation_Results']) ~= 7
    mkdir([pp,'Continuation_Results'])
end

save([pp,'Continuation_Results/contKdSteady_r',num2str(r),'v',num2str(v),'d',num2str(round(1000*D))],'kdGrid','stateKdContGrid');

%% Plotting
%{
cols = parula([length(initState)]);                              %Colors for plot

figure('visible','off')

for ii = 1 :length(initState)
    txt = ['Y', num2str(ii-1)];
    semilogx(kdGrid,20000*stateKdContGrid(ii,:),'linewidth',1,'color',cols(ii,:),'DisplayName',txt)
    hold on
end
grid on
legend show
axis([-inf inf 0 20000*M])
xlabel('K_D')
ylabel('st. state')
title(['r = ',num2str(r),', v = ',num2str(v),', d = ',num2str(D)])
pause(.01)
hold off

%}
end