% Variable latente z_n {speakers}
% Variable observada x_n {diccionario}

N = 25;	 % Numero de 'palabras' en diccionario
K = 5;	 % Numero de speakers

T = 100;  % Numero de muestras en el tiempo
EX = 2; % Numero de ejemplos

NK = N/K; % Numero de palabras por speaker
TK = T/K; % Numero de muestras por speaker


%%%%%%%%%%%% Generar cadena sintética %%%%%%%%%%%%
% Un renglon para cada speaker
data = [];
dd = round(unifrnd(1, K, K, TK)) + repmat(0:K:TK, TK, 1)';
dd = reshape(dd', T, 1)';
data = [data; dd];

% dd = round(unifrnd(1, K, K, TK)) + repmat(0:K:TK, TK, 1)';
% dd = reshape(dd', T, 1)';
% data = [data; dd];


%%%%%%%%%%%% Parámetros originales %%%%%%%%%%%%%%
orig.priori = zeros(N, 1);
orig.priori(data(1)) = 1; 
orig.priori = normalize(orig.priori, 1);

orig.mtrans = zeros(N, N);
for i = 1:N
	ii = find(data(1:T-1) == i);
	for j = 1:N-1
		orig.mtrans(i, j) = sum(data(ii+1) == j);
	end	
end
orig.mtrans = normalize(orig.mtrans, 2);

orig.memisn = zeros(K, N);
for i = 0:K:(N-1)
	ii = i + (1:NK);
	jj = (i/K) + 1;
	orig.memisn(ii, jj) = 1 / NK;
end
orig.memisn = normalize(orig.memisn, 1);


%%%%%%%%%%%% Parámetros iniciales (estim) %%%%%%%%%%%%%%
esti.priori = normalize(rand(N, 1));
esti.mtrans = normalize(rand(N, N), 2);
esti.memisn = normalize(rand(N, K), 2);

[LL, priori, mtrans, memisn, iter] = hmm_em(data, orig.priori, esti.mtrans, esti.memisn, 5)
