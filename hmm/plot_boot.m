%%%%%%%%%%%%%%%%%%%%%%
resol = '-r400';
for index = 1:(length(seq_boot))
    set(0, 'DefaultAxesFontSize', 13)
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

    
    sp(2) = plot([zvalue; zvalue], [min(yy); 1.1*max(yy)], '--r', 'LineWidth', 2);          
    
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);    
    
    %sum(dd > zvalue), ...
    tt = title(sprintf('Prueba de hipótesis para m_%d vs m_%d interlocutores',...
                seq_offs+seq_boot(index), ...
                seq_offs+seq_boot(index)+1));
    box off;
    
    set(tt, 'FontSize', 16)
    
    % print('-dpng', sprintf(''), resol);
end