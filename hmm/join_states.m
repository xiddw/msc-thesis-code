function [ param ] = join_states(param, i, j)
    n = size(param.mtrans, 1);
    
    if i > n || j > n
        return
    end
    
    idx = (param.hid == i);
    param.hid(idx) = j;
    
    mm = max(param.hid);
    
    if mm == n
        oo = setdiff(1:6, unique(param.hid));
        
        idx = (param.hid > oo);
        
        param.hid(idx) = param.hid(idx) - 1;
    end
end

