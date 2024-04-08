function y = ifngProd(x,a,b,c)
    x=x/1000;
    if nargin == 1
        a = 8.05995678728349;
        b = 0.388299803866324;
        c = 5.93007440145389;
    elseif nargin == 2
        %a = 8.05995678728349;
        b = 0.388299803866324;
        c = 5.93007440145389;
    end
    
    %% Manuela's Activation Function
    logy = max(log(15),a*(1-exp(-b*(x-c))));
    y = exp(logy);%/(5e4);
    
    %% Hill Function
    %logy = max(log(15), a*(x^b/(x^b + c^b)));
    %y = exp(logy)/(5e4);
    
end