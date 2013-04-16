%{
% cd 'C:\Users\Estudiante\Dropbox\~ESCUELA\CIMAT\4 Semestre\ST2\prog\'
cd 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\'
addpath('hmm\')
addpath('voicebox\')
addpath('voice\')
mex -O -outdir hmm hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
%}

ruta = 'pruebas/prb5/';

tic;

MAX_ITER_ESTIM = 30; % 50
MAX_ITER_HMM = 200; % 300

R_SERIES = 200; % 200

% Variable latente z_n {speakers}
% Variable observada x_n {diccionario}
K = 120;	% Numero de 'palabras' en diccionario

T = 7000;  	% Numero de muestras en el tiempo
EX = 1; 	% Numero de ejemplos

NN = 6;     % Numero de speakers

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

%%% Parámetros originales (suavizados) (corresponden al modelo 1)
oorig = params_gen(N1, K, KN1);
oorig.obs = hmm_sample(oorig, N1, K, T);
oorig.hid = key1(oorig.obs);

orig = oorig;
orig.hid = ceil(medfilt1(orig.hid, 5));

listLL1 = [];
listLL2 = [];
listbic = [];
listfp1 = [];
listfp2 = [];

for qqq = 1:10

NN = 2 + qqq;% + qqq;     % Numero de speakers

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
tic;

maxLL1 = -1e12;
maxLL2 = -1e12;

maxi1 = 0;
maxi2 = 0;

% Iterar 50 veces para cada modelo, estimando parámetros
% y conservar parámetros que correspondan a mayor verosimilitud.

fprintf('Initial parameter estimation: \n');

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

archivo = strcat(ruta, 'lists.mat');
save(archivo, 'listLL1', 'listLL2', 'listbic', 'listfp1', 'listfp2')

end

return;


% Simular bootstrapped series con modelo parametrizado (modelo 1)
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
    
    maxLL1 = -1e12;
    maxLL2 = -1e12;

    maxi1 = 0;
    maxi2 = 0;
    
    tic;    
    for ii = 1:MAX_ITER_ESTIM
        esti1 = params_rnd(N1, K, KN1);
        esti1.obs = hmm_sample(esti1, N1, K, T);
        esti1.hid = key1(esti1.obs); 

        [LL1, param, iter1] = ...
            hmm_em(boot.obs, esti1.priori, esti1.mtrans, esti1.memisn, MAX_ITER_HMM);

        if(LL1(end) > maxLL1)
            maxLL1 = LL1(end);
            maxi1 = ii;

            fin1 = param;
        end     

        esti2 = params_rnd(N2, K, KN2);
        esti2.obs = hmm_sample(esti2, N2, K, T);
        esti2.hid = key2(esti2.obs);

        [LL2, param, iter2] = ...
            hmm_em(boot.obs, esti2.priori, esti2.mtrans, esti2.memisn, MAX_ITER_HMM);

        if(LL2(end) > maxLL2)
            maxLL2 = LL2(end);
            maxi2 = ii;

            fin2 = param;
        end 

        fprintf('Iter: %2d; ', ii);
        fprintf('c: %f, m1 (%02d): %f; ', LL1(end), maxi1, maxLL1);
        fprintf('c: %f, m2 (%02d): %f; ', LL2(end), maxi2, maxLL2);
        fprintf('\n');      
    end   
    
    llrb = maxLL2 - maxLL1;
    fprintf('log LR (boot){%3d}: [%f] \n', i, llrb);
    
    if llrb > llro
        b = b+1;
    end   
    
    if maxLLRB < llrb
        maxLLRB = llrb;
    end
       
    toc;
end

pvalue = (b + 1) / (R_SERIES + 1);
