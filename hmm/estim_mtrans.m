function [ mtrans ] = estim_mtrans( grd )
    nsp = max(grd);
    lgd = length(grd);
    spk = zeros(nsp, nsp);
    
    for i = 1:(lgd-1)
        n = grd(i);
        m = grd(i+1);
        
        spk(n, m) = spk(n, m)+1;
    end
    
    mtrans = normalize(spk, 2);
end

