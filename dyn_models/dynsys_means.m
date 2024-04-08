function m = dynsys_means(kf,kb,v,vh,Tm,pia)

    % This function returns at steady state the mean value for each of the
    % variables in a contact area system, as well as the mean residency
    % time of such a NP.
    % Arguments
    % - kf(double) : forward rate of the system
    % - kb(double) : backward rate of the system
    % - v(integr) : valence of the NP considered
    % - Tm         : maximum number of contact area TCR that we want to
    % consider. (e.g. for TM=20, will compute steady state values for all
    % possible # of TCR in contact area from 1 to Tm)

    m = zeros(Tm,3);
    m(:,1) = (1:Tm)';
    m(1,2) = 1;%kf*v/(kb+kf*v);
    m(1,3) = 1/kb;
    for i = 1:Tm
        temp    = dynsys_contact(kf,kb,v,vh,i,pia);
        if ~isempty(temp.ss) && ~isempty(temp.rt)
            j = min(i,vh);
            m(i,1)  = j;                                            %Max # of pMHC-TCR complexes
            distbound = temp.ss(2:end)/sum(temp.ss(2:end));         %Distribution of bound TCR
            m(i,2)  = (1:j)*distbound;                              %Average number of bound TCR (for max i # TCR on CA)
            m(i,3)  = [temp.rt]'*distbound;                         %Average dwell time (for max i # TCR on CA)
        else
            m(i,2) = 1;%kf*v/(kb+kf*v);
            m(i,3) = 1/kb;
        end
    end
end