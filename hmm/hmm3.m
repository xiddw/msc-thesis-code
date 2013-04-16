%{

cd 'E:\ESCUELA\CIMAT\3° Semestre\PT1\data\'
addpath('hmm\')
addpath('voice\')
%kc = csvread('grab22.csv');
clc
%}

tic;

% Variable latente z_n {speakers}
% Variable observada x_n {diccionario}
N = 5;	 	% Numero de speakers
K = 100;	% Numero de 'palabras' en diccionario

T = 5000;  	% Numero de muestras en el tiempo
EX = 1; 	% Numero de ejemplos

KN = K/N; 	% Numero de palabras por speaker
TN = T/N; 	% Numero de muestras por speaker

key = reshape(repmat(1:N, KN, 1), 1, K);

%%%%%%%%%%%% Parametros originales %%%%%%%%%%%%%%
orig.priori = zeros(N, 1);
orig.priori(1) = 1; 
orig.priori = normalize(orig.priori, 1);

a = 1;
b = 200;
c = 2;
orig.mtrans = toeplitz([b, a, zeros(1, N-2)], ...
                       [b, c, zeros(1, N-2)]);                   
orig.mtrans(1, N) = a;
orig.mtrans(N, 1) = c;
           
orig.mtrans = normalize(orig.mtrans, 2);

orig.memisn = zeros(N, K);
for i = 1:N 
    orig.memisn(i, (i-1)*KN + (1:KN)) = 100;    
end

orig.memisn = orig.memisn + abs(ceil(5*randn(N, K)));
orig.memisn = normalize(orig.memisn, 2);

%%%%%%%%%%%% Generar cadena aleatoria  (orig) %%%%%%%%%%%%
data = zeros(1, T);
spkr = randsample(N, 1, true, orig.priori);
for i = 1:T
    word = randsample(K, 1, true, orig.memisn(spkr, :));
    spkr = randsample(N, 1, true, orig.mtrans(spkr, :));
    data(i) = word;
end

orig.data = key(data);


%%%%%%%%%%%% Generar cadena aleatoria  (esti) %%%%%%%%%%%%
% Un renglon para cada speaker
%{
data = [];
dd = round(unifrnd(1, KN, N, TN)) + repmat(0:KN:(K-1), TN, 1)';
dd = reshape(dd', T, 1)';
data = [data; dd];
%}
%%%%%%%%%%%% Parámetros iniciales (estim) %%%%%%%%%%%%%%
esti.priori = normalize(rand(N, 1));
esti.mtrans = normalize(rand(N, N), 2);
esti.memisn = normalize(rand(N, K), 2);
esti.data = key(data);

[LL, priori, mtrans, memisn, gamma, iter] = ...
    hmm_em(data, esti.priori, esti.mtrans, esti.memisn, 350);

fin.priori = priori;
fin.mtrans = mtrans;
fin.memisn = memisn;

gg = repmat(max(gamma), N, 1);
fin.data = mod(find(gamma == gg), N) + 1;
fin.data = fin.data';

toc;

myplot(orig, esti, fin);
%myplot(orig, fin);

print -dpng p1.png