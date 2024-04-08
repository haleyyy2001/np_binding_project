% Contact Area Dynamics Plots

pp = ['louis_results/'];
valence20 = [9,13,61,210];
valence14 = [8,11,14,31];
valence_all = [valence14; valence20];
k=1;

for r=[14,20]
    
    valence = valence_all(k,:);

    kd_range = [10:25];

    load([pp,'insprobNP_r',num2str(r),'.mat'])

    maxTCRcov = 1;
    for i = 1:length(distCov)
        maxTCRcov = max(maxTCRcov,length(distCov{i}));
    end

    % 2. Compute the unbinding time and the average number of bound TCR per NP
    %    avgs_contAreaModel = dynsys_means(koff/kd,koff,v,vh,maxTCRcov,insProbTCRtoNP);
    %   backwardRate = unbindingRate(distCov,avgs_contAreaModel(:,3));

    DwellTime = [];

    for v = valence

        load([pp,'insprob_r',num2str(r),'_v',num2str(v),'.mat'])

        vh = round(v/7);
        aveT = [];
        stStates = [];
        for kd = kd_range
            kon=1;
            koff = kd;

            for i=1:maxTCRcov
                temp = dynsys_contact(kon,koff,v,vh,i,insProbTCRtoNP);
            end

            stStates = [stStates, temp.ss];
            aveT = [aveT, mean(temp.rt)];
        end

        % Plot contact area steady states
        f1 = figure(1);
        subplot(2,4,find(valence==v)+4*(k-1))
        hold on
        for i=1:length(stStates(:,1))
            plot(kd_range, stStates(i,:), 'LineWidth',2, 'DisplayName',['X',num2str(i-1)])
        end
        hold off
        axis([10 25 0 1.0])
        xlabel('K_D')
        title(['r=',num2str(r),' v=',num2str(v)])

        set(gca,'FontSize',16,'LineWidth',2)

        DwellTime = [DwellTime; aveT];

    end

    if k==2
        leg=legend('Position', [0.035 0.5 0.01 0.1]);
        title(leg,'Bound TCRs')
    
        ax1 = subplot(2,4,1);
        ax2 = subplot(2,4,5);

        ylabel(ax1,'Steady States')
        ylabel(ax2,'Steady States')

        AddLetters2Plots({ax1,ax2},{'A','B'},'Hshift',-0.08,'Vshift',-0.07,'FontSize',22)

        f1.Units='centimeters';
        f1.Position=[1 10 50 30];

    end
    
    f2=figure(2);
    subplot(1,2,k)
    hold on
    for i = 1:length(valence)
        plot(kd_range, DwellTime(i,:),'LineWidth',2,'DisplayName',['v=',num2str(valence(i))])
    end
    hold off
    xlabel('K_D')
    title(['r=',num2str(r),' nm'])
    
    leg=legend('Location','northeast');
    title(leg, 'Valence')
    
    k=k+1;
    
end

ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

ylabel(ax1,'Average Dwell Time')

set(ax1,'Fontsize',20)
set(ax2,'Fontsize',20)

AddLetters2Plots({ax1,ax2},{'','B'},'Hshift',-0.08,'Vshift',-0.07,'FontSize',32)
AddLetters2Plots({ax1,ax2},{'A',''},'Hshift',-0.1,'Vshift',-0.07,'FontSize',32)

f2.Units='centimeters';
f2.Position=[1 50 30 15];