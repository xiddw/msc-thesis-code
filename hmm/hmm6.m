%{
    cd 'C:\Users\xiddw\Documents\GitHub\msc-thesis-code\'
    % cd 'C:\Users\Estudiante\Documents\GitHub\msc-thesis-code\'
    addpath('hmm\')
    % addpath('voicebox\')
    addpath('mfcc\')
    addpath('voice\')
mex -O -outdir hmm hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
%}

kk = [45:15:90, 100:20:200];

%grnd = 'mfcc\calderon5_ground.csv';

grnd = 'pruebas\cuervo1f_ground.csv';
%grnd = 'mfcc\noct1f_ground.csv';

kk = [120];

for www = kk
    tic;

    MAX_ITER_ESTIM = 30;
    MAX_ITER_HMM = 340;

    R_SERIES = 200;

    % Variable latente z_n {speakers}
    % Variable observada x_n {diccionario}
    
    %noct1f
    
    ruta = strcat('pruebas\prb_cuervo2f_', int2str(kk), '\')
    arch = strcat('pruebas\cuervo1f_', int2str(kk), '.csv')
    
    %ruta = strcat('mfcc\prb_noct1f_', int2str(kk), '\')
    %arch = strcat('mfcc\noct1f_', int2str(kk), '.csv')
       
    disp(strcat('Iniciando bootstrap para ->   ', arch));
    
    mkdir(ruta);

    kc = csvread(arch);

    listLL1 = [];
    listLL2 = [];
    listbic = [];
    listfp1 = [];
    listfp2 = [];

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
    
    seq_boot = 1:8;

    for qqq = seq_boot
        NN = 1 + qqq;
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

        lambda = 1.0;
        bic = maxLL1 - 0.5 * lambda * (N1-1)+(N1*(N1-1))+(N1*(K-1)) * log(T) ;

        ffin1 = sort_params(orig, fin1);
        ffin2 = sort_params(orig, fin2);

        toc; 
        img1 = strcat(ruta, int2str(NN), 'to', int2str(NN+1));
        img2 = strcat(ruta, int2str(NN), 'to', int2str(NN+1), '_');
        [fp1, fp2] = myplot(orig, img1, ffin1, img2, ffin2);

        listLL1 = [listLL1, maxLL1];
        listLL2 = [listLL2, maxLL2];

        listbic = [listbic, bic];
        listfp1 = [listfp1, fp1];
        listfp2 = [listfp2, fp2];

        archivo = strcat(ruta, int2str(NN), 'to', int2str(NN+1));
        save(archivo, 'orig', 'fin1', 'fin2', 'ffin1', 'ffin2')

    end

    archivo = strcat(ruta, 'lists.mat');
    save(archivo, 'listLL1', 'listLL2', 'listbic', 'listfp1', 'listfp2')
    
    close all;

end

N = 120;
% T = 7219;
T = 6415;

ii = 2;
hh = length(listLL1);
bb = zeros(1, hh);
MM = 1;
lambda = 3e3;
for i = 1:hh
    MM = MM+1;
    bb(i) = (listLL1(i) - listLL2(i)) - 0.5 * lambda * (MM-1)+(MM*(MM-1))+(MM*(N-1)) * log(T);
end

f0 = figure; 

sp(1) = plot(ii:(hh+1), bb, 'b');
hold on;
sp(2) = plot(ii:(hh+1), bb, 'or', 'MarkerFaceColor', 'r');
t1 = title('Selección de modelo con BIC');

set(t1, 'FontSize', 16)
set(gca, 'FontSize', 12);
set(gca, 'box', 'off');
set(sp(2),'MarkerSize', 7);
set(sp, 'linewidth', 2)

legend(strcat('lambda= ', int2str(lambda)))
hold off;