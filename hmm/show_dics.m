function show_dics(fdic, fgrd, ppath)
    if nargin < 3
        ppath = 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\';
    end
    
    ext = '.csv';
    % kk = [10:10:50, 100:100:500];
    kk = [10:10:50, 100:20:200];
    %kk = [100];
    
    for k = 1:numel(kk)
        sdic = strcat(ppath, fdic, int2str(kk(k)), ext);
        dic = csvread(sdic);
        
        sgrd = strcat(ppath, fgrd, ext);
        grd = csvread(sgrd);
        
        memisn = estim_memisn(grd, dic);        
        
        figure;
        imagesc(memisn);
    end
end


