function [LL, nmtrans, nmpriori_ini, nmemisn, nmpriori_fin, gamma] = ...
	calc_values(priori, mtrans, memisn, data)
	
	[N, K] = size(memisn);
	[EX, T] = size(data);	
	
	nmtrans = zeros(N, N);
	nmpriori_ini = zeros(N, 1);
	nmpriori_fin = zeros(N, 1);
	nmemisn = zeros(N, K);
	
	LL = 0;
		
	for ii = 1:EX
		zz = data(ii, :);		
		
		%{
		ems_like = zeros(N, T);        
		for tt = 1:T
			ems_like(:, tt) = memisn(:, zz(tt));
		end
		%}

		%%% Fordward-Backward
		%[alpha, beta, gamma, llc, xi] = fwd_bwd(priori, mtrans, ems_like);
        [~, ~, gamma, xi, llc] = cfwd_bwd(priori, mtrans, memisn, zz);

		LL = LL + llc;
		
		nmtrans = nmtrans + xi;
		
		nmpriori_ini = nmpriori_ini + gamma(:, 1);
		nmpriori_fin = nmpriori_fin + gamma(:, T);
		
		if T < K
			for t = 1:T
				k = zz(t);
				nmemisn(:, k) = nmemisn(:, k) + gamma(:, t);
			end
		else
			for k = 1:K				
                jj = (zz == k);				
                if any(jj)  
                    nmemisn(:, k) = nmemisn(:, k) + sum(gamma(:, jj), 2);
                end
            end
        end        
    end
    
% jj = find(zz == k);
% if ~isempty(jj)
		
		
		
		
