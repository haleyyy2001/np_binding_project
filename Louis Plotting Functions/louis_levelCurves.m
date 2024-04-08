% Level curves of IFNg activation ploted against NP dose and NP valence

VALENCE = [8,11,14,31,54];
DOSE = [2:30];
%VALENCE = [54,31,14,11,8];

pp = ['louis_results/'];                                % Must import IPs for all valences considered.

%data = [repmat(14,28,1), get_ifndata(14)];
data = [[repmat(20,21,1), get_ifndata(20)]];

chain = load([pp,'../fig_nsim15000date20-Jan-2023/chain1.txt']);
sschain = load([pp,'../fig_nsim15000date20-Jan-2023/sschain1.txt']);

%load('./temp_MCMC/chain_81472.mat');
%load('./temp_MCMC/sschain_81472.mat');

params = chain(find(sschain==min(sschain)),:);

if length(params) == 5
    kd = params(1); %24.7252;
    mu = params(2); %8.7209;
    alph = params(3); %0.3883;
    C = params(4); %4.3013;
    KD = params(5); %0.6;
    
    pt=1;
    
    par_lbls = {'kd','mu','alph','C','KD'};
else
    kd = params(1); %24.7252;
    pt = params(2); %0.5;
    mu = params(3); %8.7209;
    alph = params(4); %0.3883;
    C = params(5); %4.3013;
    KD = params(6); %0.6;
    
    par_lbls = {'kd','pt','mu','alph','C','KD'};
end

ifng_params = [mu, alph, C];

dtcr = 10;
ntcr = 20;
rnbh = 10;
r_ind=1;

for r=[14,20]
    y = zeros(length(VALENCE),length(DOSE));

    for i = [1:length(VALENCE)]
        v = VALENCE(i);

        if exist([pp,'insprob_r',num2str(r),'_v',num2str(v),'.mat'])==0
            louis_comp_functions(r,v,45,dtcr,rnbh,ntcr,kd);
        end

        for j = [1:length(DOSE)]
            d = DOSE(j);

            np_params = [r,v,d];

            if exist([pp,'Continuation_Results/contKdSteady_r',num2str(r),'v',num2str(v),'d',num2str(d*1000),'.mat'])==0
                louis_continuation(r,v,d,45,pp);
            end

            [y(i,j),m,mu] = louis_activationFunction(kd,KD,45,np_params,ifng_params,pt,pp);
        end
    end

    fig = figure(2);
    subplot(2,2,r_ind);
    contourf(DOSE,VALENCE,y,'ShowText','on','LineWidth',2)
    %yline(8,'r','Linewidth',1.5);
    %yline(11,'r','Linewidth',1.5);
    %xline(20,'b','Linewidth',1.5);
    %xline(10,'b','Linewidth',1.5)
    title(['($r=',num2str(r),', d_{\rm TCR}=',num2str(dtcr),'$)'],'Interpreter','Latex')
    grid on
    
    subplot(2,2,r_ind+2);
    if r==14
        v1=7;
        v2=10;
        v3=4;
        v4=13;
    else
        v1=8;
        v2=12;
        v3=4;
        v4=16;
    end
    
    hold on    
    plot(DOSE,y(v1,:),'LineWidth',2,'DisplayName',['v=',num2str(VALENCE(v1))]);
    plot(DOSE,y(v2,:),'LineWidth',2,'DisplayName',['v=',num2str(VALENCE(v2))]);
    %plot(DOSE,y(v3,:),'LineWidth',2,'DisplayName',['v=',num2str(VALENCE(v3))]);
    %plot(DOSE,y(v4,:),'LineWidth',2,'DisplayName',['v=',num2str(VALENCE(v4))]);
    hold off
    legend('Location','SouthEast')
    
    r_ind = r_ind+1;
end

ax1 = subplot(2,2,1);
yline(8,'--m','Linewidth',1.5);
yline(11,'--r','Linewidth',1.5);
ylabel('Valence','Interpreter','Latex')

ax2 = subplot(2,2,2);
yline(9,'--m','Linewidth',1.5);
yline(13,'--r','Linewidth',1.5);


ax3 = subplot(2,2,3);
xlabel('NP Concentration ($ \times 10^{11}$ NP/ml)','Interpreter','Latex')
ylabel('IFN$\gamma$ (ng/ml)','Interpreter','Latex')
grid on

ax4 = subplot(2,2,4);
xlabel('NP Concentration ($ \times 10^{11}$ NP/ml)','Interpreter','Latex');
grid on

AddLetters2Plots({ax1,ax2,ax3,ax4},{'A','B','C','D'},'Hshift',-0.08,'Vshift',-0.07,'FontSize',22) 

cb = colorbar('Location','EastOutside','FontSize',16);
title(cb,'IFN$\gamma$ (ng/ml)','Interpreter','Latex')
set(ax1,'FontSize',16)
set(ax2,'Fontsize',16)
set(ax3,'FontSize',16)
set(ax4,'Fontsize',16)
fig.Units='centimeters';
fig.Position=[5 1 30 10];
saveas(fig,['louis_results/dTCR10/Level_Curves/levelcurves_r',num2str(r),'.jpg'])