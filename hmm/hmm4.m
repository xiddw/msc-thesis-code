% Variable latente z_n {speakers}
% Variable observada x_n {diccionario}

%{
cd 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\'
addpath('hmm\')
addpath('voicebox\')
addpath('voice\')
mex -O -outdir hmm hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp

mex -g hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
kc = csvread('grab22.csv');
%}

MAX_ITER_HMM = 600;

kc = csvread('mfcc\calderon3_100.csv');

N = 4;          % Numero de speakers
K = max(kc);    % Numero de 'palabras' en diccionario

if max(size(kc)) == size(kc, 2)
    kc = kc';
end

T = size(kc, 1);  	% Numero de muestras en el tiempo

if mod(T, 2) == 1
    T = T - 1;
    kc = kc(1:T);
end

EX = 1; 	% Numero de ejemplos

KN = int32(K/N); 	% Numero de palabras por speaker
TN = int32(T/N); 	% Numero de muestras por speaker

%%%%%%%%%%%% Generar cadena sintetica %%%%%%%%%%%%
data = kc';

%%%%%%%%%%%% Parametros originales %%%%%%%%%%%%%%
orig = params_gen(N, K, KN);

orig.obs = data;
orig.hid = csvread('mfcc\calderon3_ground.csv');

orig.priori = estim_priori(orig.hid);
orig.memisn = estim_memisn(orig.hid, orig.obs);

tic;
%%%%%%%%%%%% Parametros iniciales (estim) %%%%%%%%%%%%%%
esti.priori = normalize(rand(N, 1));
esti.mtrans = normalize(rand(N, N), 2);
esti.memisn = normalize(rand(N, K), 2);

[LL, fin, iter] = ...
    hmm_em(data, esti.priori, esti.mtrans, esti.memisn, MAX_ITER_HMM);

%f = fin;

toc;
myplot(orig, 'plot1.png', fin, 'plot2.png');
