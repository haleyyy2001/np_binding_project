function [parameterGrid, stateGrid] = continuation(parameterIndex,parameterMinimum,...
    parameterMaximum,parameterStep,konInitial,koffInitial,stateInitial,tolerance,...
    D,v,vh,insProbTCRtoNP,distCov,Mvec)
    % ---------------------------------------------------------------------
    % This function runs the continuation for the steady states, for a
    % given range of parameters and a function
    %  - konInitial:       initial value of the parameter kon
    %  - koffInitial:      initial value of the parameter koff
    %  - parameterMaximum: upper bound for continuation
    %  - parameterMinimum: lower bound for continuation
    %  - parameterStep:    exponential step
    %  - stateInitial:     initial state of the system (startpoint Newton)
    %  - parameterIndex:   index of the parameter in question
    %           - 1 -> wrt kon
    %           - 2 -> wrt koff
    %           - 3 -> wrt kd = koff/kon
    %  - tolerance:        tolerance for newton Method
    %  - 
    % ---------------------------------------------------------------------
    if parameterIndex == 1
        parameterInitial = konInitial;
    elseif parameterIndex == 2
        parameterInitial = koffInitial;
        R = D*v*konInitial;
    elseif parameterIndex == 3
        parameterInitial = koffInitial/konInitial;
        R = D*v*konInitial;
    end
    
    if parameterInitial > parameterMaximum || parameterInitial < ...
            parameterMinimum
        error('Initial parameter is not within specified bounds')
    end
    
    if parameterInitial/parameterStep < parameterMinimum
    % If the parameter can not get any smaller, only move forward (j=2)
        jopts = 2;
    elseif parameterInitial*parameterStep > parameterMaximum
    % If the parameter can not get any larger, only move backward (j=1)
        jopts=1;
    else jopts = 1:2;
    end
    
    stateGrid       = stateInitial(:);
    parameterGrid   = parameterInitial;
    
    for jop = 1:length(jopts)
        j = jopts(jop);
        
        
        % Determine initial parameter in the continuation
        if j == 1
            parameterNow = parameterInitial/parameterStep;
        elseif j == 2
            parameterNow = parameterInitial*parameterStep;
        else
            error('What?')
        end
        statePrevious = stateInitial(:)./sum(stateInitial(:));
        
        
        
        
        while parameterNow < parameterMaximum && parameterNow > parameterMinimum
            if parameterIndex == 1
                R    = D*parameterNow*v;
                temp = dynsys_means(parameterNow,koffInitial,v,vh,10,insProbTCRtoNP);
                B    = unbindingRate(distCov,temp(:,3));
            elseif parameterIndex == 2
                temp = dynsys_means(konInitial,parameterNow,v,vh,10,insProbTCRtoNP);
                B    = unbindingRate(distCov,temp(:,3));
            elseif parameterIndex == 3
                temp = dynsys_means(konInitial,konInitial*parameterNow,v,vh,40,insProbTCRtoNP);
                B    = unbindingRate(distCov,temp(:,3));
            else error('Continuation can not be done - No parameter correspondence')
            end
            dX = 1      ;
            while norm(dX) > tolerance
                dX = newtonZeros(R,Mvec,B,statePrevious);
                stateNow = (statePrevious+dX)/sum(statePrevious+dX);
                statePrevious = stateNow;
            end
            
            if j == 1
                stateGrid       = [stateNow(:),  stateGrid];
                parameterGrid   = [parameterNow, parameterGrid];
                parameterNow    = parameterNow/parameterStep;
            elseif j == 2
                stateGrid       = [stateGrid,     stateNow(:)];
                parameterGrid   = [parameterGrid, parameterNow];
                parameterNow    = parameterNow*parameterStep;
            end
        end
        
        
    end
    
    
    
end