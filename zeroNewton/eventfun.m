function [x,isterm,dir] = eventfun(R,Mvec,B,y,t)
    dy = f(R,Mvec,B,t,y);
    x = norm(dy) - 1e-3;
    isterm = 1;
    dir=0;
end