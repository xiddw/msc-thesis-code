function [fp1, fp2] = myplot(objects, archivo)
    cm = colormap(hsv(10));
    
    colormap('default');
    set(0, 'DefaultAxesFontSize', 28);
    tit_fs = 30;
    
    resol = '-r400';
    
    titles = fieldnames(objects);
    nrows = 7;
    ncols = numel(titles);
        
	sf = figure;
    sp = zeros(ncols*nrows, 1);
    tt = zeros(ncols);
    
    k = 1;  
    for i = 1:ncols
        j = 1;
        name = titles{i};
        
        if i > 1
            o = objects.(titles{1});
            t = objects.(name);
            t = sort_params(o, t);
        else
            t = objects.(name);            
        end
        
        mrows = nrows;
        %sp(k) = subplot(mrows, );                        
        figure;
        stem(t.priori, 'filled', 'LineWidth', 4, 'MarkerSize', 6);
        mm = numel(t.priori);
        box off;        
        set(gca,'XTick',1:mm)
        set(gca,'XLim',[0.5, mm+0.5])
        set(get(gcf,'CurrentAxes'), 'LineWidth', 4);
        
        xlabel('Interlocutores', 'FontSize', tit_fs);
        
        if i == 1 
            ylabel('Prob. a priori', 'FontSize', tit_fs);
        end
        
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [8 4]);
        set(gcf, 'PaperPosition', [0 0 8 4]);
        
        print('-dpng', sprintf('%s%c_%d_%d', archivo, 'p', i, j), '-r400');
        
        %continue;
        
        j = j+1;
        k = k+1;

        %sp(k) = subplot(mrows, ncols, (j-1)*ncols + i);        
        figure;
        imagesc(t.mtrans);      
        
        [spk, ~] = size(t.mtrans);        
        set(gca,'XTick', 1:spk);
        set(gca,'YTick', 1:spk);
        xlabel('Interlocutores', 'FontSize', tit_fs);
        ylabel('Interlocutores', 'FontSize', tit_fs);

        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [8 5]);
        set(gcf, 'PaperPosition', [0 0 8 5]);
        
        print('-dpng', sprintf('%s%c_%d_%d', archivo, 'p', i, j), '-r400');

        sy = max(t.memisn(:));
        sy = str2double(sprintf('%5.2f', sy)) + 0.01;
        s = 1;
        
        %cm = colormap(jet(spk));
        
        for s = 1:spk        
            [sa, sb] = size(t.memisn);
            % sp(k) = subplot(mrows, ncols, (j-1)*ncols + i);
            figure;    
            
            stem(t.memisn(s, :), 'filled', 'Color', cm(s, :), 'LineWidth', 4, 'MarkerSize', 5);
            box off;
            set(gca,'XLim', [1, sb]);
            set(gca,'YLim', [0, sy]);
            set(get(gcf,'CurrentAxes'), 'LineWidth', 4);
            xlabel('Diccionario', 'FontSize', tit_fs);
            ylabel('Probabilidad', 'FontSize', tit_fs);        

            j = j+1;
            k = k+1;
            
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [8 3]);
            set(gcf, 'PaperPosition', [0 0 8 3]);
            
            print('-dpng', sprintf('%s%c_%d_%d', archivo, 'p', i, j), '-r400');
        end      
    end
    
    % Graficas de cambio entre speakers

    %% 
    close all;
    return;
    
    nrows = ncols;
    ncols = 1;   
    
    grd = titles{1};
    a = objects.(grd);
    
    for j = 2:nrows
        sp = [];
        tt = [];
    
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
        mm = max(a.hid);
        ll = length(a.hid);
        set(gca,'YTick', 1:mm);
        set(gca,'YLim', [0.5, mm+0.5]);
        set(gca,'XLim', [1, ll]);
        la(1) = xlabel('Tiempo (ms)');
        la(2) = ylabel('Interlocutores');
        
        i = i+1;    

        sp(i) = subplot(2, ncols, i);
        hold on;        
        sl(i) = plot(b.hid, '-b');
        mm = max(b.hid);
        ll = length(b.hid);
        set(gca,'YTick', 1:mm);
        set(gca,'YLim', [0.5, mm+0.5]);
        set(gca,'XLim', [1, ll]);
        la(3) = xlabel('Tiempo (ms)');
        la(4) = ylabel('Interlocutores');
        
        p = abs(a.hid - b.hid);
        p(p == 0) = nan;    
        plot(p, '.r');

        fp1 = sum(p > 0);    
        fprintf('(Modelo %d) FP + FN = %d [%0.5f]\n', ...
            max(b.hid), fp1, fp1/length(a.hid));

        i = i+1;

        set(sp, 'FontSize', 13)
        set(sp, 'box', 'off')
        set(sp, 'color', 'white')
        set(sp, 'linewidth', 1)        
        set(la, 'FontSize', 15)

        set(sl, 'linewidth', 2)
        set(tt, 'FontSize', 16)           
        
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [16 5]);
        set(gcf, 'PaperPosition', [0 0 16 5]);  
               
        print('-dpng', sprintf('%s%c_%d', archivo, 's', j), '-r400');    
    end
    %%%
    
    close all;

    %%%

    
    %% set(sf, 'Position', get(0,'Screensize')); 
    %% saveas(sf, archivo2);    

    
    %colormap(pink)