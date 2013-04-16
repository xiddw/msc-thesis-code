function R = close_to_zero(a, b, tol)
    
    if nargin < 2
        b = 0;
    end
    
    if nargin < 3 
        tol = 1e-6;
    end

	R = sum(abs(a - b) < tol);    