function [R, f] = normalize(A, dim)
    if nargin < 2
        f = sum(sum(A, 1));
        f = f + (f == 0);
        R = A ./ f;		
	else
		f = sum(A, dim);		
		f = f + (f == 0);
		
		if dim == 1
			dims = [size(A, 1), 1];
		else
			dims = [1, size(A, 2)];
		end
		
		ff = repmat(f, dims);
		R = A ./ ff;
    end
end 

% f = sum(A(:));
% f = f + (f == 0);
% R = A ./ f;