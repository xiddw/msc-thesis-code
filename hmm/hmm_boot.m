%{
    cd 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\'
    
    cd 'C:\Users\xiddw\Documents\GitHub\msc-thesis-code\'
    %cd 'C:\Users\Estudiante\Documents\GitHub\msc-thesis-code\'
    addpath('hmm\')
    addpath('mfcc\')
    addpath('voice\')
mex -O -outdir hmm hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
%}

kk = [45:15:90, 100:20:200];

stem = 'cuervo1f';

grnd = strcat('pruebas\', stem, '_ground.csv');

MAX_ITER_ESTIM = 30;
MAX_ITER_HMM = 340;

R_SERIES = 700;

kk = 140;

seq_boot = 3:6;
seq_offs = 1;
ss = length(seq_boot);

for www = kk
    tic;
    % Variable latente z_n {speakers}
    % Variable observada x_n {diccionario}
       
    ruta = strcat('pruebas\prb_b2_', stem, '_', int2str(kk), '\');
    arch = strcat('pruebas\', stem, '_', int2str(kk), '.csv');
    disp(ruta);disp(arch);
    
    %ruta = strcat('mfcc\prb_noct1f_', int2str(kk), '\')
    %arch = strcat('mfcc\noct1f_', int2str(kk), '.csv')
       
    disp(strcat('Iniciando bootstrap para ->   ', arch));
    
    mkdir(ruta);

    kc = csvread(arch);

    K = max(kc);    % Numero de 'palabras' en diccionario
    % NN = 2;          % Numero de speakers

    if max(size(kc)) == size(kc, 2)
        kc = kc';
    end

    T = numel(kc);  	% Numero de muestras en el tiempo

    if mod(T, 2) == 1
        T = T - 1;
        kc = kc(1:T);
    end

    data = kc';
    
    listLL1 = zeros(1, ss);
    listLL2 = zeros(1, ss);
    listLLR = zeros(1, ss);
    listpva = zeros(1, ss);
    %listbic = zeros(1, ss);
    listfp1 = zeros(1, ss);
    listfp2 = zeros(1, ss);
    
    listLLRB = zeros(ss, R_SERIES);

    for qqq = seq_boot
        NN = seq_offs + qqq;
        
        %%% Primer modelo
        N1 = NN;	% Numero de speakers
        KN1 = int32(K/N1); % Numero de palabras por speaker
        TN1 = int32(T/N1); % Numero de muestras por speaker

        key1 = reshape(repmat(1:N1, KN1, 1), 1, KN1 * N1);
        key1 = [key1, repmat(key1(end), 1, K - numel(key1))];

        %%% Segundo modelo
        N2 = N1+1;	% Numero de speakers + 1
        KN2 = int32(K/N2); % Numero de palabras por speaker
        TN2 = int32(T/N2); % Numero de muestras por speaker

        key2 = reshape(repmat(1:N2, KN2, 1), 1, KN2 * N2);
        key2 = [key2, repmat(key2(end), 1, K - numel(key2))];

        %%% Parámetros originales (corresponden al modelo 1)
        orig = params_gen(N1, K, KN1);
        orig.obs = data;
        orig.hid = csvread(grnd);

        orig.priori = estim_priori(orig.hid);
        orig.memisn = estim_memisn(orig.hid, orig.obs);

        tic;

        maxLL1 = -1e12;
        maxLL2 = -1e12;

        maxi1 = 0;
        maxi2 = 0;

        % Iterar MAX_ITER_ESTIM veces para cada modelo, estimando parámetros
        % y conservar parámetros que correspondan a mayor verosimilitud.

        for ii = 1:MAX_ITER_ESTIM
            esti1 = params_rnd(N1, K, KN1);
            esti1.obs = hmm_sample(esti1, N1, K, T);
            esti1.hid = key1(esti1.obs); 

            [LL1, param iter1] = ...
                hmm_em(orig.obs, esti1.priori, esti1.mtrans, esti1.memisn, MAX_ITER_HMM);

            if(LL1(end) > maxLL1)
                maxLL1 = LL1(end);
                maxi1 = ii;        
                fin1 = param;
            end      

            esti2 = params_rnd(N2, K, KN2);
            esti2.obs = hmm_sample(esti2, N2, K, T);
            esti2.hid = key2(esti2.obs);

            [LL2, param, iter2] = ...   
                hmm_em(orig.obs, esti2.priori, esti2.mtrans, esti2.memisn, MAX_ITER_HMM);

            if(LL2(end) > maxLL2)
                maxLL2 = LL2(end);
                maxi2 = ii;        
                fin2 = param;
            end 
            
            %fprintf('Iter: %2d; ', ii);
            %fprintf('c1: %f, m1 (%02d): %f; ', LL1(end), maxi1, maxLL1);
            %fprintf('c2: %f, m2 (%02d): %f; ', LL2(end), maxi2, maxLL2);
            %fprintf('\n');   
            fprintf('------------\n');
        end

        % Estimar log LikelihoodRatio (Observed)
        llro = maxLL2 - maxLL1;
        fprintf('maxLL1: [%f], maxLL2: [%f]\n', maxLL1, maxLL2);
        fprintf('log LR (obs): [%f] \n', llro);

        toc; 
        
        %%
        ffin1 = sort_params(orig, fin1);
        ffin2 = sort_params(orig, fin2);
        
        
        %%
        fprintf('Bootstrapped series: \n');

        b = 0;

        maxLLRB = -1e12;

        for i = 1:R_SERIES    
            fprintf('\tSeries %3d: \n', i);

            %{
            fin1.obs = hmm_sample(fin1, N1, K, T);
            [ll1] = calc_values(fin1.priori, fin1.mtrans, fin1.memisn, fin1.obs);

            fin2.obs = hmm_sample(fin2, N2, K, T);      
            [ll2] = calc_values(fin2.priori, fin2.mtrans, fin2.memisn, fin2.obs);
            %}

            boot = ffin1;
            boot.obs = hmm_sample(boot, N1, K, T);

            tmaxLL1 = -1e12;
            tmaxLL2 = -1e12;

            maxi1 = 0;
            maxi2 = 0;

            tic;    
            for ii = 1:(MAX_ITER_ESTIM)
                esti1 = params_rnd(N1, K, KN1);
                esti1.obs = hmm_sample(esti1, N1, K, T);
                esti1.hid = key1(esti1.obs); 

                [LL1, param, iter1] = ...
                    hmm_em(boot.obs, esti1.priori, esti1.mtrans, esti1.memisn, MAX_ITER_HMM);

                if(LL1(end) > tmaxLL1)
                    tmaxLL1 = LL1(end);
                    maxi1 = ii;

                    fin1 = param;
                end     

                esti2 = params_rnd(N2, K, KN2);
                esti2.obs = hmm_sample(esti2, N2, K, T);
                esti2.hid = key2(esti2.obs);

                [LL2, param, iter2] = ...
                    hmm_em(boot.obs, esti2.priori, esti2.mtrans, esti2.memisn, MAX_ITER_HMM - 50);

                if(LL2(end) > tmaxLL2)
                    tmaxLL2 = LL2(end);
                    maxi2 = ii;

                    fin2 = param;
                end 
                
                fprintf('*');
                % fprintf('\t\tIter: %2d; ', ii);
                % fprintf('c: %f, m1 (%02d): %f; ', LL1(end), maxi1, tmaxLL1);
                % fprintf('c: %f, m2 (%02d): %f; ', LL2(end), maxi2, tmaxLL2);
                % fprintf('\n');      
            end   

            llrb = tmaxLL2 - tmaxLL1;
            fprintf('\nraxLL1: [%f], raxLL2: [%f]\n', tmaxLL1, tmaxLL2);
            fprintf('\tlog LR (boot){%3d}: [%f] \n', i, llrb);
            
            listLLRB(qqq, i) = llrb;

            if llrb > llro
                b = b+1;
            end   

            if maxLLRB < llrb
                maxLLRB = llrb;
            end

            toc;
        end

        pvalue = (b + 1) / (R_SERIES + 1);

        %%
        %%{       
        img1 = strcat(ruta, int2str(NN), 'to', int2str(NN+1));
        img2 = strcat(ruta, int2str(NN), 'to', int2str(NN+1), '_');
        [fp1, fp2] = myplot(orig, img1, ffin1, img2, ffin2);
        %%}
        
        %%
        listLL1(qqq) = maxLL1;
        listLL2(qqq) = maxLL2;
        listLLR(qqq) = maxLLRB;
        % listbic(qqq) = bic;
        listfp1(qqq) = fp1;
        listfp2(qqq) = fp2;
        
        listpva(qqq) = pvalue;
        
        archivo = strcat(ruta, int2str(NN), 'to', int2str(NN+1));
        save(archivo, 'orig', 'fin1', 'fin2', 'ffin1', 'ffin2')
        
        %% FOR THE LOLZ (en caso de que se interrumpa ejecución o algo u,u)
        N = max(orig.hid);
        archivo = strcat(ruta, 'lists.mat');
        save(archivo, 'K', 'T', 'N', 'seq_offs', 'seq_boot', ...
                  'listLL1', 'listLL2', 'listLLR', 'listLLRB',  ...
                  'listfp1', 'listfp2', 'listpva');        

    end

    N = max(orig.hid);
    
    archivo = strcat(ruta, 'lists.mat');
    save(archivo, 'K', 'T', 'N', 'seq_offs', 'seq_boot', ...
                  'listLL1', 'listLL2', 'listLLR', 'listLLRB',  ...
                  'listfp1', 'listfp2', 'listpva');
    
    close all;
    toc;
    
    %%%%%%%%%%%%%$$$$$$$$$$$$$$
    %%%%%%%%%%%%%$$$$$$$$$$$$$$
    resol = '-r400';
    
    K = max(orig.obs);
    N = max(orig.hid);
    T = length(orig.obs);
    
    lll = (1:1e2:50e2);
    hh = length(listLL1);
    bb = zeros(hh, length(lll));
    j = 1;
    for l = lll
        MM = seq_offs;
        lambda = l;
        for i = 1:hh
            MM = MM+1;
            bb(i, j) = (listLL1(i)) - 0.5 * lambda * (MM-1)+(MM*(MM-1))+(MM*(K-1)) * log(T);
        end
        %figure; plot(seq_offs+(1:hh), bb(:, j));
        j = j+1;
    end
    figure; 
    set(0, 'DefaultAxesFontSize', 13)
    
    surfc(lll, seq_offs+(1:hh), bb); colormap('cool');
    tt = title(sprintf('Superficie de curvas BIC para distintos valores de lambda'));
    set(tt, 'FontSize', 16)
    
    x = 1:length(lll);
    xx = lll;
    y = seq_offs+(1:hh);
    z = bb;
    
    figure;
    [px, py] = gradient(z, 1, 2);
    pc = contour(x, y, z, 30, 'LineWidth', 1); 
    colormap('cool');
    hold on; 
    quiver(x, y, px, py, 'k', 'LineWidth', 1);
    
    st = max(diff(xx));
    %xt = str2num(get(gca, 'XTickLabel'));
    %set(gca,'XTickLabel', xt*st);
    
    tt = title(sprintf('Curva de nivel de superficie BIC para distintos valores de lambda'));
    set(tt, 'FontSize', 16)
    
end

        