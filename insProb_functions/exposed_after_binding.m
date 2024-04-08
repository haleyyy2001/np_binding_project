function y = exposed_after_binding(pos,sam,dis)
    % Variables:
    % - pos: positions of points,
    % - sam: sampled points,
    % - dist: distance from sampled points.
    % This function returns the positions of all points in [pos] that are at
    % a distance greater than [dis] from the points in [sam].

    
    d = dist(sam',pos);
    min_d = min(d,[],1);
    
    y = pos(:,find(min_d > dis));
end