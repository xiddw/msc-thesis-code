function [ param ] = sort_params(orig, dest, collapse)    
    if nargin < 3
        collapse = false;
    end
    
    if numel(orig.hid) < numel(dest.hid)
        [orig.hid, dest.hid] = resize_vectors(orig.hid, dest.hid);
    else
        [dest.hid, orig.hid] = resize_vectors(dest.hid, orig.hid);
    end
    
    offset = 100;    
    n = numel(unique(dest.hid));
    
    difsize = (numel(unique(orig.hid)) ~= numel(unique(dest.hid)));
    
    per = zeros(n, 1);
    % Ordering hidden states.
    if(max(orig.hid) ~= max(dest.hid))
        for i = 1:n
            per(i) = mode(orig.hid(dest.hid == i));
        end
    else
        for i = 1:n
            per(i) = mode(orig.hid(dest.hid == i));
            % [~, per] = max(dest.mtrans, [], 1);
        end      
    end
    
    %% Code to collapse when an item has more states than other one
    if difsize && collapse
        dd = dist(per');
        dd = dd + diag(nan(1, n));
        
        unique(dest.hid);

        [x, y] = find(dd == min(dd(:)), 1, 'first');

        dest.hid(dest.hid == y) = x;
        
        u = unique(dest.hid);
        
        if max(u) == n
            dest.hid(dest.hid == n) = (min(u) - 1);
        end    
        
        [param] = sort_params(orig, dest);
    else
    % normal code
        [simi, pos] = pdist2(orig.memisn, dest.memisn, ...
                                'cosine', 'smallest', 1);
                        
        [~, per] = sort(pos + simi);
        %[~, per] = sort(per);

        param = dest;
        
        param.memisn = param.memisn(per, :);
        
        for i = 1:n
            idx = (param.hid == i);
            param.hid(idx) = per(i) + offset;
        end

        param.hid = param.hid - offset;
        
        [~, per] = max(param.mtrans, [], 1);
        param.mtrans = param.mtrans(per, :);   
        param.priori = param.priori(per);   
    end
    
    
end

