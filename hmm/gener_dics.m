function gener_dics( )
    
    %kk = [10:10:50, 100:20:200];
    kk = [45:15:90, 100:20:200];
    
    %numc = 5;
    
    fpath = 'E:/ESCUELA/CIMAT/4 Semestre/ST2/prog/mfcc/';
    ffile = 'cuervo1f';
    fext  = '.csv';
    
    disp(strcat('Input: ', ffile, fext));
    
    for k = 1:length(kk)
       
        disp(strcat('Generando diccionario para k = ', int2str(kk(k))));
        
        finput = strcat(fpath, ffile, '_ceps', fext);
        foutput = strcat(fpath, ffile, '_', int2str(kk(k)), fext);
        
        aa = csvread(finput)'; %%%%%%%%%
        
        %{
        bb = aa - mean(aa, 1)
        [W, ~] = princomp(bb);
        size(W);
        %}
        
        bb = fkmeans(aa, kk(k));        
        csvwrite(foutput, bb);        
    end
end

