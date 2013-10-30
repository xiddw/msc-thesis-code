%%%%%%%%%%%%%%%%%%%%%%
resol = '-r100';
typef = '-depsc';

font = 16;
tit_fs = 18;
for index = seq_boot
    set(0, 'DefaultAxesFontSize', font)
    dd = listLLRB(index, :);
    nd = length(dd);
    [yy, xx] = ksdensity(dd, 'npoints', nd);
    figure; 
    
    zvalue = listLL2(index) - listLL1(index);
    ii = int16(0.95*nd);
    signif = xx(ii);    
    
    xmin = 0.9*min(xx);
    xmax = 1.1*max([xx, zvalue]);
    
    ymin = min(yy);
    ymax = 1.1*max(yy);    
    
    rectangle('Position', [signif, ymin, xmax, ymax], ...
              'FaceColor', [225, 240, 225]/255, ...
              'LineStyle', 'none');              
          
    hold on; 
        
    plot(xx, yy, 'LineWidth', 2);
    xlabel('log-likelihood ratio observada', 'FontSize', tit_fs);
    ylabel('Densidad de prob.', 'FontSize', tit_fs);
    
    plot([zvalue; zvalue], [min(yy); 1.1*max(yy)], '--r', 'LineWidth', 3);
    
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);    
    
    box off;
    
    comp = '<';
    if zvalue > signif; comp = '>'; end;
    
    fprintf('(Modelo %d) zvalue : %9.5f %c [reject: %9.5f]\n', ...
            index, zvalue, comp, signif);    
    
    if exist('archivo', 'var') 
        name = archivo;
        if exist('farchivo', 'var') 
            name = farchivo;
        end
        
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [9 3]);
        set(gcf, 'PaperPosition', [0 0 9 3]);
        
        arch = strcat(name, 'boot', int2str(index));
        print(typef, arch, resol); 
        close;
    end    
end