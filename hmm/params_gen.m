function [ param ] = params_gen( N, K, KN )
    param.priori = zeros(N, 1);
    param.priori(1) = 1; 
    param.priori = normalize(param.priori, 1);

    a = 1;
    b = 100;
    c = 5;
    param.mtrans = toeplitz([b, a, zeros(1, N-2)], ...
                           [b, c, zeros(1, N-2)]);                   
    param.mtrans(1, N) = a;
    param.mtrans(N, 1) = c;

    param.mtrans = normalize(param.mtrans, 2);

    param.memisn = zeros(N, K);
    
    for i = 1:N
        param.memisn(i, (i-1)*KN + (1:KN)) = 100;    
    end
    
    param.memisn = param.memisn(1:N, 1:K);
    
    param.memisn = param.memisn + abs(ceil(5*randn(N, K)));
    param.memisn = normalize(param.memisn, 2);
end

