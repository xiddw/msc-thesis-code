function [ va, vb, vc ] = resize_vectors( va, vb, vc )
% Ajusta el tamaño de dos vectores para que concuerden
% el vector va lo hace del tamño del vector vb 
% (para esto numel(va) < numel(vb))
    if nargin < 3
        vc = vb;
    end
        
    nr = round(numel(vb) / numel(va));
    
    aa = repmat(va, [nr 1]);
    
    va = aa(:)';
    
    la = numel(va);
    lb = numel(vb);
    
    if la < lb
        vb = vb(1:la);
        vc = vc(1:la);
    elseif la > lb
        va = va(1:lb);
    end
end

    