function [ data ] = hmm_sample( param, N, K, T )
    data = zeros(1, T);
    
    % rand sample vector (N x 1) with replacement and weighted
    spkr = randsample(N, 1, true, param.priori);
    
    
    for i = 1:T
        word = randsample(K, 1, true, param.memisn(spkr, :));
        spkr = randsample(N, 1, true, param.mtrans(spkr, :));
        data(i) = word;
    end

    % params.data = key(data);
end