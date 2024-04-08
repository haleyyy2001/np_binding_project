function m = dynsys_means_v2_markovchain(kf,kb,vi,Tm,pia)

    % This function returns at steady state the mean value for each of the
    % variables in a contact area system, as well as the mean residency
    % time of such a NP.
    % Arguments
    % - kf(double) : on rate of the system
    % - kb(double) : off rate of the system
    % - vi(integr) : valence of the NP considered
    % - Tm         : maximum number of contact area TCR that we want to
    % consider. (e.g. for TM=20, will compute steady state values for all
    % possible # of TCR in contact area from 1 to Tm)
    
    m = zeros(Tm,3);
    m(:,1) = (1:Tm)';
    m(1,2) = 1;%kf*vi/(kb+kf*vi);
    m(1,3) = 1/kb;
    for i = 1:Tm
        temp   = dynsys_contact_v2_markovchain(kf,kb,vi,i,pia);
        if ~isempty(temp.ss) && ~ isempty(temp.rt)
            j = min(i,length(temp.ss));
            m(i,1)  = j;
            distengage = temp.ss;
            m(i,2)  = (1:j)*distengage;
            m(i,3)  = temp.rt;
        end
    end
end