function [param] = params_rnd(N, K, ignor)
    param.priori = normalize(rand(N, 1));
    param.mtrans = normalize(rand(N, N), 2);
    param.memisn = normalize(rand(N, K), 2);
end

