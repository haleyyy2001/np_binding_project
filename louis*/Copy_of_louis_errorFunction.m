function er = louis_errorFunction(params,data,ansr)
%   Error function for a set of parameters

    a_comp = zeros(length(data(:,1)),1);
    a_call = zeros(length(data(:,1)),1);
    mAvg = zeros(length(data(:,1)),1);
    mDis = zeros(length(data(:,1)),1);
    er = 0;  % Included to initialize er for individual fitting.
    
    % Allow for choice of parameters to test, depending on the length of
    % the parameter vector
    
    kd          = params(1);
    %pho_params  = params(2);       %[params, 5]; %params(6:end);
    pho_params = 1;               % When reparametrizing with A = pt*alph and C* = C/pt
    ifng_params = params(2:4);
    %ifng_params = [params(3),1,1];
    KD          = params(5);
    
    
    theta_deg   = 45;               %round(params(5))*5 + 20; % 7 angles to choose from from 20 to 50
    
    pp = ['louis_results/'];
    
       
    for i = 1:length(data(:,1))
        r = data(i,1);
        v = data(i,2);
        d = data(i,3);
        dind = i;
        a_d = data(i,4);
        
        np_params = [r,v,d];
                  
        [a_call(i), mAvg(i), mDis(i)] = louis_activationFunction(kd,KD,theta_deg,np_params,ifng_params, pho_params,pp);
        
        if i > 1
            j=0;
            while ((a_call(i) < a_call(i-1)) && (v == data(i-1,2)) && (j < 5))
                [a_call(i), mAvg(i), mDis(i)] = louis_activationFunction(kd,KD,theta_deg,np_params,ifng_params,pho_params,pp);
                j = j+1;
            end
        end
        
        
        a_comp(i) = (a_d-a_call(i))^2;
        
        % For fitting curves individually
        %if (r == 14) && (v == 8)
        %    er = er + a_comp(i);
        %end
         
    end
    er = (sum(a_comp));        % For fitting all curves together
    
    if nargin == 3  
        if ansr == 'y';
            data = [data,a_call,mAvg,mDis];
            rs = unique(data(:,1));
            for i = 1:length(rs)
                datar = data(data(:,1)==rs(i),:);
                vs = unique(datar(:,2), 'first');
    
                cols = cool(round(1.3*length(vs)));
                
                hold off
                for j=1:length(vs)
                    datav = datar(datar(:,2)==vs(j),:);
                    
                    disp(datav(:,4));
                    
                    txt1 = ['Exp - v',num2str(vs(j))];
                    txt2 = ['Mod - v',num2str(vs(j))];
                 
                    figure(i)
                    semilogx(datav(:,3),datav(:,4),':o','color',cols(j,:),'linewidth',2, 'DisplayName', txt1)
                    hold on
                    semilogx(datav(:,3),datav(:,5),'-s','color',0.7*cols(j,:),'linewidth',2, 'DisplayName', txt2)
                    %{
                    figure(i)
                    semilogx(datav(:,3),datav(:,6),'-o','color',cols(j,:),'linewidth',1.5)
                    hold on
                    semilogx(datav(:,3),datav(:,7),'--o','color',cols(j,:),'linewidth',1.5)
                    %}
                end
                
                figure(i)
                title(['Activation of r=', num2str(rs(i)),' NPs'])
                xlabel('NP dose (\times 10^{11})')
                ylabel('IFN_\gamma (ng/mL)')
                legend('Location', 'northeastoutside')
                legend show
                grid on
                set(gca,'FontSize',16,'LineWidth',2)
                
                %figure(i)
                %legend('averages','gamma')
            end
        end
    end
end
