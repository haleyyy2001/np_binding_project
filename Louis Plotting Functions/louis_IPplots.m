% Insertion probabilities plots for CAM and TNM.
% Also including distribution of covered TCRs per NP.
%% Contact Area Insertion Probabilities

pp = ['louis_results/'];%/dTCR10/New IP/

valence14 = [8, 11, 14, 31, 54];
valence20 = [9, 13, 61, 210];
valence = {valence14, valence20};

r = 20;
r_ind = 1;
for r = [14,20];
    
    f1 = figure(1);
    subplot(1,2,r_ind);
    hold on

    cols = 0.8*parula(length(valence{r_ind}));
    v_i = 1;
    for v = valence{r_ind}
        load([pp,'insprob_r',num2str(r),'_v',num2str(v),'.mat'])
        plot([0:length(insProbTCRtoNP)],[insProbTCRtoNP,0],'color',cols(v_i,:),'DisplayName',[num2str(v)],'LineWidth',2)
        v_i=v_i+1;
    end

    hold off
    xlabel('Bound pMHC-TCR complexes')
    title(['r=',num2str(r),' nm'])
    leg = legend();
    title(leg, 'Valence')
    set(gca, 'FontSize', 16, 'LineWidth', 2)
    
    r_ind=r_ind + 1;
end

ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

ylabel(ax1,'Insertion Probability')
xlim(ax1,[0 3]);
xticks(ax1, [0 1 2 3])
xlim(ax2,[0 6]);

AddLetters2Plots({ax1,ax2}, {'A',''},'Hshift',-0.06,'Vshift',-0.075,'FontSize',22)
AddLetters2Plots({ax1,ax2},{'','B'},'Hshift',-0.06,'Vshift',-0.075,'FontSize',22)

f1.Units='centimeters';
f1.Position=[1 10 35 10];



    %% Nanocluster Insertion Probabilities
    
r_ind = 1;
    
for r=[14,20];
    
    load([pp,'insprobNP_r',num2str(r),'.mat'])

    f2 = figure(2);
    subplot(1,2,r_ind);
    hold on
    yyaxis left
    plot([0:length(insProbNPtoNC)], [insProbNPtoNC,0],'LineWidth',2)
    ylabel(['Insertion Probability'],'Interpreter','Latex')

    txt = '$max(c_{NP})$';

    yyaxis right
    plot([1:length(probBindNPi)], probBindNPi,'LineWidth',2)
    xline(length(probBindNPi),'--','Linewidth',2)
    ylabel('TCR$_{\rm nc}$ with $c_{NP} \geq i$','Interpreter','Latex')
    text(length(probBindNPi)-3.8,0.8,txt,'Fontsize',16,'Interpreter','Latex');
    axis([0 20 0 1])
    xlabel('Number of NPs bound to Nanocluster','Interpreter','Latex')
    title(['r=',num2str(r),' nm'])
    set(gca,'FontSize',22,'LineWidth',2)
    hold off

    r_ind = r_ind +1;
end

ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

AddLetters2Plots(f2, {'A',''},'Hshift',-0.105,'Vshift',-0.075,'FontSize',32)
AddLetters2Plots({ax1,ax2},{'','B'},'Hshift',-0.06,'Vshift',-0.045,'FontSize',32)
%annotation(f2,'textbox',[0 0.1 0 0.1],'String','A','EdgeColor','None','Fontsize',22)
%annotation(f2,'textbox',[0.1 0.2 0.7 0.8],'String','B','EdgeColor','None','Fontsize',32)

f2.Units='centimeters';
f2.Position=[1 10 50 15];


%% Nanocluster Model - Distribution of covered TCRs
r_i =1;
for r=[20]

    load([pp,'insprobNP_r',num2str(r),'.mat'])

    mx = 1;
    f3=figure(3);
    %subplot(1,2,r_i)
    cols = parula(length(distCov)/2);
    hold on
    for i=0:floor(length(distCov)/2)-1
        temp = length(distCov{1+2*i});
        mx = max(mx, temp);

        plot([1:temp],distCov{1+2*i},'DisplayName',num2str(1+2*i),'Linewidth',2)

    end

    hold off
    leg = legend('Location','Northeast');
    title(leg,'Order of NP','Fontsize',12)
    title(['r=',num2str(r),' nm'])
    xlabel('TCRs covered by NP')
    
    set(gca,'FontSize',20,'Linewidth',2','colororder',cols)

    r_i = r_i+1;
end

f3.Units='centimeters';
f3.Position=[1 10 45 15];

%ax1 = subplot(1,2,1);
%ax2 = subplot(1,2,2);
xlim([1 3])
%xlim(ax2,[1 4])

ylabel('Proportion of Simulated TCR_{nc}');
%{
AddLetters2Plots({ax1,ax2},{'A',''},'Hshift',-0.07,'Vshift',-0.075,'FontSize',32) 
AddLetters2Plots({ax1,ax2},{'','B'},'Hshift',-0.07,'Vshift',-0.065,'FontSize',32) 
%}