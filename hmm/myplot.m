function [fp1, fp2] = myplot(objects, archivo)
    cm = ametrine(7);
    
    colormap('default');
    set(0, 'DefaultAxesFontSize', 28);
    tit_fs = 30;
    
    resol = '-r100';
    typef = '-depsc';
    
    titles = fieldnames(objects);
    nrows = 7;
    ncols = numel(titles);
        
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
        stem(t.priori, 'filled', 'Marker', 's', 'LineWidth', 6, 'MarkerSize', 8);
        mm = numel(t.priori);
        box off;        
        set(gca,'XTick',1:mm)
        set(gca,'XLim',[0.5, mm+0.5])
        set(get(gcf,'CurrentAxes'), 'LineWidth', 3);
        
        xlabel('Interlocutores', 'FontSize', tit_fs);
        ylabel('Prob. a priori', 'FontSize', tit_fs);
        
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [8 4]);
        set(gcf, 'PaperPosition', [0 0 8 4]);
        
        print(typef, sprintf('%s%c_%d_%d', archivo, 'p', i, j), resol);
        
        %continue;
        
        j = j+1;
        k = k+1;

        %sp(k) = subplot(mrows, ncols, (j-1)*ncols + i);        
        figure;
        imagesc(t.mtrans);      
        wp = 8;
        if i == ncols
            colorbar();
            wp = 9;
        end
        colormap(ametrine);
        
        [spk, ~] = size(t.mtrans);        
        set(gca,'XTick', 1:spk);
        set(gca,'YTick', 1:spk);
        xlabel('Interlocutores', 'FontSize', tit_fs);
        ylabel('Interlocutores', 'FontSize', tit_fs);

        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [wp 5]);
        set(gcf, 'PaperPosition', [0 0 wp 5]);
        
        print(typef, sprintf('%s%c_%d_%d', archivo, 'p', i, j), resol);

        sy = max(t.memisn(:));
        sy = str2double(sprintf('%5.2f', sy)) + 0.01;
        s = 1;
        
        %cm = colormap(jet(spk));
        
        for s = 1:spk        
            [sa, sb] = size(t.memisn);
            % sp(k) = subplot(mrows, ncols, (j-1)*ncols + i);
            figure;    
            
            stem(t.memisn(s, :), 'filled', 'Marker', 's', 'Color', cm(s, :), 'LineWidth', 6, 'MarkerSize', 8);
            box off;
            set(gca,'XLim', [1, sb]);
            set(gca,'YLim', [0, sy]);
            set(get(gcf,'CurrentAxes'), 'LineWidth', 3);
            xlabel('Diccionario', 'FontSize', tit_fs);
            ylabel('Probabilidad', 'FontSize', tit_fs);        

            j = j+1;
            k = k+1;
            
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [8 3]);
            set(gcf, 'PaperPosition', [0 0 8 3]);
            
            print(typef, sprintf('%s%c_%d_%d', archivo, 'p', i, j), resol);
        end      
    end
    
    % Graficas de cambio entre speakers

    %% 
    
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
        
        %%%%%%%%%%%% STARTS SEQ %%%%%%%%%%%%

        %sp(i) = subplot(2, ncols, i);    
        plot(a.hid, '-b', 'LineWidth', 2);
        mm = max(a.hid);
        ll = length(a.hid);
        
        set(gca,'FontSize', 16);
        set(gca,'YTick', 1:mm);
        set(gca,'YLim', [0.5, mm+0.5]);
        set(gca,'XLim', [1, ll]);
        box off;
        
        xlabel('Tiempo (ms)', 'FontSize', 18);
        ylabel('Interlocutores', 'FontSize', 18);
        
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [16 3]);
        set(gcf, 'PaperPosition', [0 0 16 3]);  
               
        print(typef, sprintf('%s%c_%d_1', archivo, 's', j), resol);            
        
        %sp(i) = subplot(2, ncols, i);
        plot(b.hid, '-b', 'LineWidth', 2);
        mm = max(b.hid);
        ll = length(b.hid);
        
        set(gca,'FontSize', 16);
        set(gca,'YTick', 1:mm);
        set(gca,'YLim', [0.5, mm+0.5]);
        set(gca,'XLim', [1, ll]);
        box off;
        
        xlabel('Tiempo (ms)', 'FontSize', 18);
        ylabel('Interlocutores', 'FontSize', 18);
        
        hold on;
        p = abs(a.hid - b.hid);
        q = b.hid;
        q(p == 0) = nan;    
        plot(q, 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

        fp1 = sum(p > 0);    
        fprintf('(Modelo %d) FP + FN = %d [%0.5f]\n', ...
            max(b.hid), fp1, fp1/length(a.hid));
        
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [16 3]);
        set(gcf, 'PaperPosition', [0 0 16 3]);  
               
        print(typef, sprintf('%s%c_%d_2', archivo, 's', j), resol);        
    end
    %%%
    
    close all;

    %%%

    
    %% set(sf, 'Position', get(0,'Screensize')); 
    %% saveas(sf, archivo2);    

    
    %colormap(pink)