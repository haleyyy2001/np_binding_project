function er = errorFunctionParFitting(params,data)
%   Error function for a set of parameters
    
    a_comp = zeros(length(data(:,1)),1);
    a_call = zeros(length(data(:,1)),1);
    
    % Allow for choice of parameters to test, depending on the length of
    % the parameter vector
    
    lp = length(params);
    if lp==3
        kd    = params(1);
        kp    = params(2);
        kappa = params(3);
        
        pi = (kp/(kp+1))^kappa;
        
    elseif lp==4
        kd    = params(1);
        kp    = params(2);
        kappa = params(3);
        mu    = params(4);
        
        pi = (kp/(kp+1))^kappa;
        
    elseif lp==5
        kd    = params(1);
        pi    = params(2);
        mu    = params(3);
        alph  = params(4);
        C     = params(5);
    end
    
        
    for i = 1:length(data(:,1))
        r = data(i,1);
        v = data(i,2);
        d = data(i,3);
        dind = i;
        a_d = data(i,4);
        
        if lp == 3
            a_call(i) = activationFunction(r,v,d,kd,pi);
        elseif lp==4
            a_call(i) = activationFunction(r,v,d,kd,pi,mu);
        elseif lp==5
            a_call(i) = activationFunction(r,v,d,kd,pi,mu,alph,C);
        end
        
        a_comp(i) = (a_d-a_call(i))^2;
    end
    er = sum(a_comp);
    
    data = [data,a_call];
    
    rs = unique(data(:,1));
    for i = 1:2
        datar = data(data(:,1)==rs(i),:);
        vs = unique(datar(:,2));
        
        cols = cool(2*length(vs));
        figure(i+2)
        hold off
        for j=1:length(vs)
            datav = datar(datar(:,2)==vs(j),:);
            semilogx(datav(:,3),datav(:,4),'color',cols(j,:),'linewidth',1.5)
            hold on
            semilogx(datav(:,3),datav(:,5),':','color',cols(j,:),'linewidth',1.5)
        end
    end
    
end
