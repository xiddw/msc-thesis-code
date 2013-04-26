function [fp1, fp2] = myplot(a, archivo1, b, archivo2, c)
    resol = '-r400';
    nrows = 4;
    
    ns = [max(a.hid), max(b.hid)];
    if nargin < 4
        error('Se requieren al menos cuatro parametros.')
    elseif nargin == 4
        ncols = 2;        
    else
        ns(3) = max(c.hid);
        ncols = 3;
    end
    
    i = 1;
    
	sf = figure;
	sp(i) = subplot(nrows, ncols, i);    
	imagesc(a.priori)
    tt(1) = title('Ground truth');
    ll(1) = ylabel('Priori');
    i = i+1;
    
	sp(i) = subplot(nrows, ncols, i);    
	imagesc(b.priori)
    tt(2) = title(strcat('Modelo para k=', int2str(ns(2))));
    i = i+1;
    
    if ncols == 3
        sp(i) = subplot(nrows, ncols, i);                
        imagesc(c.priori)
        tt(3) = title(strcat('Modelo para k=', int2str(ns(3))));
        i = i+1;
    end

	sp(i) = subplot(nrows, ncols, i);    
	imagesc(a.mtrans)
    ll(2) = ylabel('M. Transición');
    i = i+1;    

	sp(i) = subplot(nrows, ncols, i);    
	imagesc(b.mtrans)
    i = i+1;
    
    if ncols == 3
        sp(i) = subplot(nrows, ncols, i);        
        imagesc(c.mtrans)
        i = i+1;
    end

	sp(i) = subplot(nrows, ncols, i);
    imagesc(a.memisn)
    ll(3) = ylabel('M. Emisión');
    i = i+1;

	sp(i) = subplot(nrows, ncols, i);    
	imagesc(b.memisn)
    i = i+1;
    
    if ncols == 3        
        sp(i) = subplot(nrows, ncols, i);        
        imagesc(c.memisn)
        i = i+1;
    end
    
	sp(i) = subplot(nrows, ncols, i);    
	imagesc(a.hid)
    ll(4) = ylabel('Secuencia');
    i = i+1;

	sp(i) = subplot(nrows, ncols, i);    
	imagesc(b.hid)    
    i = i+1;
    
    if ncols == 3
        sp(i) = subplot(nrows, ncols, i);
        imagesc(c.hid)
    end
    
    set(sp, 'FontSize', 13)
    set(sp, 'box', 'off')
    set(sp, 'color', 'white')
    set(sp, 'linewidth', 1)      
    
    set(tt, 'FontSize', 16)
    set(ll, 'FontSize', 15)
    
    set(sp(end-ncols:end), 'yticklabel', []);
    
    % Graficas de cambio entre speakers
    nrows = ncols;
    ncols = 1;   
    
    sp = [];
    tt = [];
    
    %% set(sf, 'Position', get(0,'Screensize')); 
    %% saveas(sf, archivo1);
    if ~strcmp(archivo1, '')
        print('-dpng', archivo1, resol);
    end
    %% 
    
    %%%
    if numel(a.hid) < numel(b.hid)
        [a.hid, b.hid] = resize_vectors(a.hid, b.hid);
    else
        [b.hid, a.hid] = resize_vectors(b.hid, a.hid);
    end
    %%%
    sf = figure;    
    
    i = 1;
    
    frac = 1;
    n = length(a.hid);
    offset = 0; int64(n / 1);
    
    idx = offset + (1:(int32(n/frac)));
    
	sp(i) = subplot(nrows, ncols, i);    
	sl(i) = plot(a.hid, '-b');
    ylim([1, ns(1)+1])
    tt(1) = title('Secuencia original');
    i = i+1;    
   
    sp(i) = subplot(nrows, ncols, i);
    hold on;        
	sl(i) = plot(b.hid, '-b');
    ylim([1, ns(2)+1])
    tt(2) = title(strcat('Secuencia recuperada para k=', int2str(ns(2))));
    
    p = abs(a.hid - b.hid);
    p(p == 0) = nan;    
    plot(p, '.r');
    
    fp1 = sum(p > 0);    
    fprintf('(Modelo 1) FP + FN = %d\n', fp1);
    
    i = i+1;
    
    if nrows == 3
        sp(i) = subplot(nrows, ncols, i);
        hold on;        
               
        sl(i) = plot(c.hid(idx), '-b')
        ylim([1, ns(3)+1])
        tt(3) = title(strcat('Secuencia recuperada para k=', int2str(ns(3))));
        
        q = abs(a.hid - c.hid);  
        q(q == 0) = nan;        
        plot(q, '.r');
        
        fp2 = sum(q > 0);        
        fprintf('(Modelo 2) FP + FN = %d\n', fp2);        
        i = i+1;
    end
    
    set(sp, 'FontSize', 13)
    set(sp, 'box', 'off')
    set(sp, 'color', 'white')
    set(sp, 'linewidth', 1)
    
    set(sl, 'linewidth', 2)
    
    set(tt, 'FontSize', 16)
    
    %% set(sf, 'Position', get(0,'Screensize')); 
    %% saveas(sf, archivo2);    
    if ~strcmp(archivo2, '')
        print('-dpng', archivo2, resol);
    end
    
    %colormap(pink)