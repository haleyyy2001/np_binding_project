function [distcov,probNPi] = distCovTCR2(r,ntcr,dtcr,simdist,mx)
    % This function simulates a large number of TCRnc to determine the
    % distribution of the number of TCR that a bound NP covers in general.
    % It returns a cell, for which element i corresponds to the
    % probability distribution of the number of covered TCR by the i-th NP
    % to be bound to the TCRnc.
    rnc = 100;
    
    covMatrix = zeros(simdist,mx);
    for ii = 1:simdist
        % Determine initial positions for TCR
        pos_tcr     = rand_scatter_disk(rnc,ntcr,dtcr); 
        i=1;            % Initialize the number of NP to bind.
        nexp = ntcr;    % Initialize the number of exposed TCR
        pos_tcrtemp = pos_tcr;
        while nexp > 0
            % Sample positions for occupied TCR
            pos_np  = sample_min_dist(pos_tcrtemp,1,2*r);
            
            % Find the ones that are in the contact area
            distances = dist(pos_np',pos_tcr);
            ncov = sum(distances<=r/sqrt(2));
            covMatrix(ii,i) = ncov;
            
            % Update the number of exposed TCR and their position
            pos_tcrtemp = exposed_after_binding(pos_tcrtemp,pos_np,2*r);  %% Should be 2*r for tcrs available for next binding (used to be r)
            if isempty(pos_tcrtemp)
                nexp =0;
            else
                nexp = length(pos_tcrtemp(1,:));
            end
            i=i+1;
        end
    end
    distcov = cell(1,mx);
    for j=1:mx
        nsucc = sum(covMatrix(:,j)>0);
        disp(nsucc);
        psucctab(j) = nsucc/simdist;
        for jj = 1:max(covMatrix(:,j))
            distcov{j}(jj) = sum(covMatrix(:,j) == jj)/nsucc;
        end
    end
    %
    probNPi = psucctab;
end