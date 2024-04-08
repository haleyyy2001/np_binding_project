function dX = newtonZeros(R,M,B,X)

    % Given R, M (1->N vector), B (1->N vector) and X (0->N vector): select all
    % column vectors
    
    N = length(M);
    
    %% Jacobian matrix of system 
    
    
    diagc = zeros(N+1,1);
    diagl = zeros(N+1,1);
    diagu = zeros(N+1,1);
    
    diagc = -R*([ (M-X(2:end))./M ; 0 ] + [0; X(1:end-1)./M] ) - [ 0 ; B ];
    diagu = [ 0 ; B ] + R*[0; X(1:end-1)./M];
    diagl = R*[ (M-X(2:end))./M ; 0 ];
    
    jacMat = spdiags([diagl,diagc,diagu],[-1,0,1],N+1,N+1);
    jacMat(end,:) = ones(1,N+1);
    
    rhs = - f(R,M,B,X);
    rhs(end)=0;
    
    %cond(jacMat)               % To check for singular matrices
    
    dX = jacMat\rhs;
    
end
