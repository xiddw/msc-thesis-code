%{
    cd 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\'
    
    cd 'C:\Users\xiddw\Documents\GitHub\msc-thesis-code\'
    addpath('hmm\')
    addpath('mfcc\')
    addpath('voice\')
mex -O -outdir hmm hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
%}

kk = [45:15:90, 100:20:200];

grnd = 'pruebas\cats1f_ground.csv';

MAX_ITER_ESTIM = 30;
MAX_ITER_HMM = 340;

R_SERIES = 200;

% kk = [160];

seq_offs = 1;
seq_boot = 1:7;    
ss = length(seq_boot);

T = 0;

for www = kk
    tic;
    % Variable latente z_n {speakers}
    % Variable observada x_n {diccionario}
       
    ruta = strcat('pruebas\prb_tt_cats1f_', int2str(www), '\')
    arch = strcat('pruebas\cats1f_', int2str(www), '.csv')
    
    disp(strcat('Iniciando BIC para ->   ', arch));
    
    mkdir(ruta);

    kc = csvread(arch);

    K = max(kc);    % Numero de 'palabras' en diccionario
    %N = max(grnd);          % Numero de speakers

    if max(size(kc)) == size(kc, 2)
        kc = kc';
    end

    T = numel(kc);  	% Numero de muestras en el tiempo

    if mod(T, 2) == 1
        T = T - 1;
        kc = kc(1:T);
    end

    data = kc';
    
    listLL1 = zeros(ss);
    listLL2 = zeros(ss);    
    listfp1 = zeros(ss);
    listfp2 = zeros(ss);
    % listbic = zeros(ss);    
    
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
            %fprintf('curr1: %f, mm1 (%02d): %f; ', LL1(end), maxi1, maxLL1);
            %fprintf('curr2: %f, mm2 (%02d): %f; ', LL2(end), maxi2, maxLL2);
            %fprintf('\n');  

        end

        % Estimar log LikelihoodRatio (Observed)
        llro = maxLL2 - maxLL1;
        fprintf('log LR (obs): [%f] \n', llro);

        % lambda = 1.0;
        % bic = maxLL1 - 0.5 * lambda * (N1-1)+(N1*(N1-1))+(N1*(K-1)) * log(T) ;

        ffin1 = sort_params(orig, fin1);
        ffin2 = sort_params(orig, fin2);

        toc; 
        img1 = strcat(ruta, int2str(NN), 'to', int2str(NN+1));
        img2 = strcat(ruta, int2str(NN), 'to', int2str(NN+1), '_');
        [fp1, fp2] = myplot(orig, img1, ffin1, img2, ffin2);

        listLL1(qqq) = maxLL1;
        listLL2(qqq) = maxLL2;

        %listbic(qqq) = bic;
        listfp1(qqq) = fp1;
        listfp2(qqq) = fp2;

        archivo = strcat(ruta, int2str(NN), 'to', int2str(NN+1));
        save(archivo, 'orig', 'fin1', 'fin2', 'ffin1', 'ffin2')

    end
    
    N = max(orig.hid);

    archivo = strcat(ruta, 'lists.mat');
    save(archivo, 'K', 'T', 'N', 'seq_offs', 'seq_boot', ...
                  'listLL1', 'listLL2', 'listfp1', 'listfp2')
    
    close all;
    toc;
    
    %%%%%%%%%%%%%$$$$$$$$$$$$$$
    %%%%%%%%%%%%%$$$$$$$$$$$$$$
    
    lll = (1:2e1:40e2);
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
    figure; surfc(lll, seq_offs+(1:hh), bb); colormap cool;
end

resol = '-r400';
lp = 100;

figure;
sp(1) = plot(seq_offs+(1:hh), bb(:, lp), '-b');
hold on;
sp(2) = plot(seq_offs+(1:hh), bb(:, lp), 'or', 'MarkerFaceColor', 'r');
tt = title(sprintf('Selección BIC con lambda=%d',lll(lp)-1));
box off;
set(sp, 'linewidth', 2)
set(0, 'DefaultAxesFontSize', 13)
set(tt, 'FontSize', 16)

print('-dpng', 'cats11.png', resol);

index = 2;
set(0, 'DefaultAxesFontSize', 13)
a = ksdensity(listLLRB(index, :));
figure; 
sp(1) = plot(a);
hold on; 
sp(2) = plot([listLLR(index); listLLR(index)], [min(a); 1.1*max(a)], '--r');
tt = title(sprintf('Prueba de hipótesis para m_%d vs m_%d interlocutores',...
            index+seq_offs+seq_boot(index), ...
            index+seq_offs+seq_boot(index+1)));
box off;
set(sp, 'linewidth', 2)
set(tt, 'FontSize', 16)