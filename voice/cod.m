tic;

sp = [];
sl = [];
st = [];

%[x fs bps] = wavread("E:/PT1/data/amorosos_fin.wav");
% [x fs bps] = wavread('track8a.wav');

arch1 = 'mfcc/nocturno1.wav';
arch2 = 'mfcc/nocturno1f.wav';

[x fs bps] = wavread(arch1);

xmin = min(x);
xmax = max(x);

dur = length(x)/fs;

ws = 0.1 * fs;
ws = ws + (mod(ws, 2) == 1);
gx = ones(ws, 1)';

xx = fftconv(abs(x), gx)';
xx = xx(ws/2:end-(ws/2));

tol_silence = quantile(xx, 0.10);

pos_silence = find(xx < tol_silence);

num_silence = length(pos_silence)-1;
cc = jet(num_silence);

tt = (1:length(x))/fs;

sp(1) = subplot(6, 1, 1);
sl(1) = plot(tt, x, 'k');

%st(1) = title( 'Señal de audio original'); 
ll(1) = ylabel('Amplitud');
%ll(4) = xlabel('Tiempo (s)');

sp(2) = subplot(6, 1, 2);
sl(2) = plot(tt, x, 'k');

hold on;
j = 1;
zz = zeros(length(x), 1);

pp = 1;
pq = 1;

c = 0;
for i = 1:num_silence
    a = pos_silence(i) + ws/2;    
    b = pos_silence(i+1) - ws/2;
    
    if (a > b || b - a < 3*ws); continue; end;
    
    z = x(a:b);
    
    pq = b - a;

    tt = (a:b) / fs;
	
    cl(1) = plot([a, a]/fs, [xmin, xmax]);
    cl(2) = plot([b, b]/fs, [xmin, xmax]);
    
    bd = plot([c, a]/fs, [0, 0], '-w', 'linewidth', 2.0);
    c = b;    
    
    %st(2) = title( 'Señal de audio procesada'); 
    ll(2) = ylabel('Amplitud');
    
    set(cl, 'color', 'r', ...
            'linewidth', 3.0, ...
            'linestyle', '--');      
        
    set(bd, 'color', 'w', ...
            'linewidth', 2.0, ...
            'linestyle', '-');        
    
    j = j + 1;
    zz(pp: (pp+pq)) = z;
    
    pp = pp + pq;
end
%ll(3) = xlabel('Tiempo (s)');

zz = zz(1:pp);

wavwrite(zz, fs, arch2);

min_x = 0;
max_x = length(zz)/fs;

max_y = max(ceil(abs(10*[min(zz), max(zz)])) / 10);
min_y = -max_y;

%{
sp(3) = subplot(3, 1, 3);
zz = [zz; zeros(length(x) - length(zz), 1)];
tt = (1:length(zz))/fs;
sl(3) = plot(tt, zz, 'k');

st(3) = title( 'Señal de audio segmentada'); 
ll(4) = ylabel('Amplitud');
%}

%set(sp, 'FontSize', 14)
set(sp, 'box', 'off')
set(sp, 'color', 'white')
set(sp, 'linewidth', 1)
set(sp, 'FontSize', 12)

set(ll, 'fontSize', 12)

%set(sp([1]), 'xticklabel', []);
set(sp, 'ycolor', 'black');

%set(st, 'FontSize', 16)

set(sl, 'linewidth', 1.1)

%set(sp, 'color', 'red')

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc = melcepst(zz);
cc = flipud(cc');
% subplot(3, 1, 2);
% imagesc(cc);


[ka kb kc] = kmeans(cc', 300);
% hold on;
% subplot(3, 1, 3);

% imagesc([min_x, max_x], [min_y/4, max_y/4], kc');
% imagesc([min_x, max_x], [min_y/4, max_y/4], kc');
% plot(tt, zz, 'k');

% imagesc([min_x, max_x], [-0.01, 0.01], kc');
% plot([min_x, max_x], [min_y, max_y], tt, zz, 'k');

% hold off;
% colormap("summer");

toc