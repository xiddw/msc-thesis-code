function [ priori ] = estim_priori( grd )
    N = max(grd);
    priori = zeros(N, 1);
    priori(grd(1)) = 1; 
    %priori = normalize(orig.priori, 1);
end

