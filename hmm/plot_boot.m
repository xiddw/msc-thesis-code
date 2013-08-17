%%%%%%%%%%%%%%%%%%%%%%
resol = '-r400';
font = 14;
for index = seq_boot
    set(0, 'DefaultAxesFontSize', font)
    dd = listLLRB(index, :);
    nd = length(dd);
    [yy, xx] = ksdensity(dd, 'npoints', nd);
    figure; 
    
    zvalue = listLL2(index) - listLL1(index);
    signif = xx(int16(0.95*nd));    
    
    xmin = 0.9*min(xx);
    xmax = 1.1*max([xx, zvalue]);
    
    ymin = min(yy);
    ymax = 1.1*max(yy);    
    
    rectangle('Position', [signif, ymin, xmax, ymax], ...
              'FaceColor', [235, 255, 235]/255, ...
              'LineStyle', 'none');    
          
          
    hold on; 
        
    sp(1) = plot(xx, yy, 'LineWidth', 2);
    xlabel('log-likelihood ratio observada');
    ylabel('Densidad de probabilidad');

    
    sp(2) = plot([zvalue; zvalue], [min(yy); 1.1*max(yy)], '--r', 'LineWidth', 2);          
    
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);    
    
    %sum(dd > zvalue), ...
    %{
    tt = title(sprintf('Prueba de hipótesis para m_%d vs m_%d interlocutores',...
                seq_offs+ index, ...%%seq_boot(index), ...
                seq_offs+ index+1));%%seq_boot(index)+1));
    %}
    box off;
    
    %set(tt, 'FontSize', font+3)
    
    if exist('archivo', 'var') 
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [8 3]);
        set(gcf, 'PaperPosition', [0 0 8 3]);
        
        arch = strcat(archivo, '_', int2str(index));
        print('-dpdf', arch, resol); close;
    end    
    
    % print('-dpng', sprintf(''), resol);
end