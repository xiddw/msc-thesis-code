%{
    cd 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\'
    
    cd 'C:\Users\xiddw\Documents\GitHub\msc-thesis-code\'
    cd 'C:\Users\Estudiante\Documents\GitHub\msc-thesis-code\'
    addpath('hmm\')
    addpath('mfcc\')
    addpath('voice\')
mex -O -outdir hmm hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
%}

kk = [45:15:90, 100:20:200];

%grnd = 'mfcc\calderon5_ground.csv';

grnd = 'pruebas\cuervo1f_ground.csv';
%grnd = 'pruebas\noct1f_ground.csv';

MAX_ITER_ESTIM = 30;
MAX_ITER_HMM = 340;

R_SERIES = 200;

kk = [120];

T = 0;

for www = kk
    tic;
    % Variable latente z_n {speakers}
    % Variable observada x_n {diccionario}
       
    ruta = strcat('pruebas\prb_cuervo1f_boot_', int2str(kk), '\')
    arch = strcat('pruebas\cuervo1f_', int2str(kk), '.csv')
    
    %ruta = strcat('mfcc\prb_noct1f_', int2str(kk), '\')
    %arch = strcat('mfcc\noct1f_', int2str(kk), '.csv')
       
    disp(strcat('Iniciando bootstrap para ->   ', arch));
    
    mkdir(ruta);

    kc = csvread(arch);

    K = max(kc);    % Numero de 'palabras' en diccionario
    NN = 2;          % Numero de speakers

    if max(size(kc)) == size(kc, 2)
        kc = kc';
    end

    T = numel(kc);  	% Numero de muestras en el tiempo

    if mod(T, 2) == 1
        T = T - 1;
        kc = kc(1:T);
    end

    data = kc';
    
    seq_boot = 1:3;    
    seq_offs = 4;
    ss = length(seq_boot);
    
    listLL1 = zeros(ss);
    listLL2 = zeros(ss);
    listLLR = zeros(ss);
    listpva = zeros(ss);
    %listbic = zeros(ss);
    listfp1 = zeros(ss);
    listfp2 = zeros(ss);

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
            
            fprintf('Iter: %2d; ', ii);
            fprintf('c1: %f, m1 (%02d): %f; ', LL1(end), maxi1, maxLL1);
            fprintf('c2: %f, m2 (%02d): %f; ', LL2(end), maxi2, maxLL2);
            fprintf('\n');   
            
        end

        % Estimar log LikelihoodRatio (Observed)
        llro = maxLL2 - maxLL1;
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

                fprintf('\t\tIter: %2d; ', ii);
                fprintf('c: %f, m1 (%02d): %f; ', LL1(end), maxi1, tmaxLL1);
                fprintf('c: %f, m2 (%02d): %f; ', LL2(end), maxi2, tmaxLL2);
                fprintf('\n');      
            end   

            llrb = tmaxLL2 - tmaxLL1;
            fprintf('\tlog LR (boot){%3d}: [%f] \n', i, llrb);

            if llrb > llro
                b = b+1
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

    end

    archivo = strcat(ruta, 'lists.mat');
    save(archivo, 'listLL1', 'listLL2', 'listbic', 'listfp1', 'listfp2')
    
    close all;
    tac;
end

        