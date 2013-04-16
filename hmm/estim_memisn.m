function [ memisn ] = estim_memisn( grd, dic )
    if numel(grd) < numel(dic)
        [grd, dic] = resize_vectors(grd, dic);
    else
        [dic, grd] = resize_vectors(dic, grd);
    end

    N = numel(unique(grd));
    K = numel(unique(dic));

    memisn = zeros(N, K);

    for i = 1:N
       idx = (grd == (i));
       memisn(i, :) = hist(dic(idx), 1:K);
    end

    memisn = normalize(memisn, 2);
end

