% [LL, priori, mtrans, memisn, iter] = hmm_em(data, orig.priori, orig.mtrans, orig.memisn, 5)

function [LL, param, iter] = ...
   hmm_em(data, priori, mtrans, memisn, max_iter, verbose)
    if nargin < 6
        verbose = false;
    end
	
	iter = 0;
	
	LL = zeros(max_iter, 1);
	ll0 = -inf; 
	ll1 = 0;
	finished = 0;
	
	[N ~] = size(memisn);
	%[EX T] = size(data);
    
    %data = num2cell(data, 2);	
    while(iter < max_iter && ~finished)		
    
		%%% EM-step
		[ll1, nmtrans, ~, nmemisn, nmpriori_ini, gamma] = ...
			calc_values(priori, mtrans, memisn, data);			        
        
        finished = em_converged(ll1, ll0);        
        if(finished) 
            break;
        end
              
		priori = normalize(nmpriori_ini, 1);
		mtrans = normalize(nmtrans, 2);
		memisn = normalize(nmemisn, 2);
        
		iter = iter + 1;
		
		LL(iter) = ll0;
		ll0 = ll1;		
        
        if verbose
            fprintf('ii %d, ll = %f\n', iter, ll1);     
        end
    end
    
    LL = LL(1:iter);
    
    param.priori = priori;
    param.mtrans = mtrans;
    param.memisn = memisn;
    
    param.obs = data;
    
    gg = repmat(max(gamma), N, 1);
    param.hid = mod(find(gamma == gg), N) + 1;
    param.hid = param.hid';
	
	
