cd 'E:\ESCUELA\CIMAT\3° Semestre\PT1\data\'

%addpath('hmm_octave\')
addpath('hmm\')
addpath('voicebox\')
addpath('voice\')

mex -g hmm\cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
mex -g hmm\cfwd_bwd_2.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp

%mex -g hmm_octave\cfwd_bwd.cpp hmm_octave\cpptipos\matriz.cpp hmm_octave\cpptipos\vector.cpp

N = 5;	 	% Numero de speakers
K = 10;	 	% Numero de 'palabras' en diccionario


T = 15;  	% Numero de muestras en el tiempo
EX = 1; 	% Numero de ejemplos

KN = K/N; 	% Numero de palabras por speaker
TN = T/N; 	% Numero de muestras por speaker


data = [];
dd = round(unifrnd(1, KN, N, TN)) + repmat(0:KN:(K-1), TN, 1)';
dd = reshape(dd', T, 1)';
data = [data; dd];

esti.priori = normalize(rand(N, 1));
esti.mtrans = normalize(rand(N, N), 2);
esti.memisn = normalize(rand(N, K), 2);

memisn = esti.memisn;
zz = data;

ems_like = zeros(N, T);
        
for tt = 1:T
	ems_like(:, tt) = memisn(:, zz(tt));
end

[alpha, beta, gamma, loglike, xi] = ...
	fwd_bwd(esti.priori, esti.mtrans, ems_like);

[alpha2, beta2, gamma2, xi2, llc2] = ...
	cfwd_bwd(esti.priori, esti.mtrans, esti.memisn, data);

[alpha3, beta3, gamma3, xi3, llc3] = ...
	cfwd_bwd_2(esti.priori, esti.mtrans, esti.memisn, data);
