function [fp1, fp2] = myplot(objects, archivo1, archivo2)
    resol = '-r400';
    
    titles = fieldnames(objects);
    nrows = 4;
    ncols = numel(titles);
        
	sf = figure;
    sp = zeros(ncols*nrows, 1);
    tt = zeros(ncols);
    
    k = 1;    
    for i = 1:ncols
        j = 1;
        name = titles{i};
        t = objects.(name);

        sp(k) = subplot(nrows, ncols, (j-1)*ncols + i);
        imagesc(t.priori)
        tt(i) = title(name);
        ll(1) = ylabel('Priori');
        j = j+1;
        k = k+1;

        sp(k) = subplot(nrows, ncols, (j-1)*ncols + i);
        imagesc(t.mtrans)
        ll(2) = ylabel('M. Transición');
        j = j+1;
        k = k+1;

        sp(k) = subplot(nrows, ncols, (j-1)*ncols + i);
        imagesc(t.memisn)
        ll(3) = ylabel('M. Emisión');
        j = j+1;
        k = k+1;

        sp(k) = subplot(nrows, ncols, (j-1)*ncols + i);
        imagesc(t.hid)
        ll(4) = ylabel('Secuencia');
        j = j+1;
        k = k+1;
    end
    
    set(sp, 'FontSize', 13)
    set(sp, 'box', 'off')
    set(sp, 'color', 'white')
    set(sp, 'linewidth', 1)      
    
    %set(tt, 'FontSize', 16)
    set(ll, 'FontSize', 15)
    
    set(sp(end-ncols:end), 'yticklabel', []);
    
    % Graficas de cambio entre speakers
    sp = [];
    tt = [];
    %% 
    
    nrows = ncols;
    ncols = 1;   
    
    grd = titles{1};
    a = objects.(grd);
    
    for j = 2:nrows
        seq = titles{j};
        b = objects.(seq);
        i = j-1;
        
        ns = [max(a.hid), max(b.hid)];
        
        if numel(a.hid) < numel(b.hid)
            [a.hid, b.hid] = resize_vectors(a.hid, b.hid);
        else
            [b.hid, a.hid] = resize_vectors(b.hid, a.hid);
        end
        
        sf = figure;    
    
        i = 1;

        frac = 1;
        n = length(a.hid);
        offset = 0; int64(n / 1);

        idx = offset + (1:(int32(n/frac)));

        sp(i) = subplot(2, ncols, i);    
        sl(i) = plot(a.hid, '-b');
        ylim([1, ns(1)+1])
        %tt(1) = title(grd);
        i = i+1;    

        sp(i) = subplot(2, ncols, i);
        hold on;        
        sl(i) = plot(b.hid, '-b');
        ylim([1, ns(2)+1])
        %tt(2) = title(strcat('Secuencia recuperada para k=', int2str(ns(2))));

        p = abs(a.hid - b.hid);
        p(p == 0) = nan;    
        plot(p, '.r');

        fp1 = sum(p > 0);    
        fprintf('(Modelo 1) FP + FN = %d\n', fp1);

        i = i+1;

        set(sp, 'FontSize', 13)
        set(sp, 'box', 'off')
        set(sp, 'color', 'white')
        set(sp, 'linewidth', 1)

        set(sl, 'linewidth', 2)

        set(tt, 'FontSize', 16)
        
        if ~strcmp(archivo2, '')
            print('-dpng', archivo2, resol);
        end
    end
    %%%

    %%%

    
    %% set(sf, 'Position', get(0,'Screensize')); 
    %% saveas(sf, archivo2);    

    
    %colormap(pink)