% Variable latente z_n {speakers}
% Variable observada x_n {diccionario}

N = 3;	 % Número de 'palabras' en diccionario
K = 5;	 % Número de speakers

T = 7; 	 % Numero de muestras en el tiempo
EX = 3;  % Numero de ejemplos

%%%%%%%%%%%% Parámetros originales %%%%%%%%%%%%%%
orig.priori = normalize(rand(N, 1));	% Prob. a priori  p(z_1) 
orig.mtrans = normalize(rand(N, N), 2); % Matriz de transicion  p(z_n | z_n-1)
orig.memisn = normalize(rand(N, K), 1); % Matriz de emisiones   p(x_n | z_n) %%% 2 %%%

%%%%%%%%%%%% Parámetros iniciales (estim) %%%%%%%%%%%%%%
esti.priori = normalize(rand(N, 1));
esti.mtrans = normalize(rand(N, N), 2);
esti.memisn = normalize(rand(N, K), 1); %%% 2 %%%

%%%%%%%%%%%% Generar datos a partir de parám. originales %%%%%%%%%%%%%%

[obs, hid] = hmm_sample(orig.priori, orig.mtrans, orig.memisn, T, EX)
data = hmm_sample(orig.priori, orig.mtrans, orig.memisn, T, EX)

zz = data(1, :)
ems_like = multinomial_prob(zz, orig.memisn)
