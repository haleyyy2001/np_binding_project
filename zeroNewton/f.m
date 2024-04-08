function y = f(R,M,B,X,t)
    % R = binding rate vector
    % M = total distribution (keep at 1 to make into distribution)
    % B = backward rates
    R=R(:);M=M(:);B=B(:);X=X(:);
    y = R*( [0; X(1:end-1).*(M-X(2:end))./M] - [X(1:end-1).*(M-X(2:end))./M;0] ) +...
        [B.*X(2:end);0] - [0; B.*X(2:end)];
                
end