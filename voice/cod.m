%{
    cd 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\'
    
    cd 'C:\Users\xiddw\Documents\GitHub\msc-thesis-code\'
    %cd 'C:\Users\Estudiante\Documents\GitHub\msc-thesis-code\'
    addpath('hmm\')
    addpath('mfcc\')
    addpath('voice\')
    addpath('freezeColors\')
mex -O -outdir hmm hmm/cfwd_bwd.cpp hmm\cpptipos\matriz.cpp hmm\cpptipos\vector.cpp
%}

resol = '-r300';
typef = '-depsc';

lbl_fs = 16;
tit_fs = 14;

tic;

sp = [];
sl = [];
st = [];

%[x fs bps] = wavread("E:/PT1/data/amorosos_fin.wav");
% [x fs bps] = wavread('track8a.wav');

arch1 = 'voice/aud/amorosos_fin.wav';
arch2 = 'voice/aud/amorosos_finf.wav';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 
plot(tt, x, 'k');
box off; 
set(gca, 'FontSize', tit_fs);
ylabel('Amplitud', 'FontSize', lbl_fs);
xlabel('Tiempo (s)', 'FontSize', lbl_fs);

% Print figure to pdf and eps files    
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 3]);
set(gcf, 'PaperPosition', [0 0 8 3]);

print('-dpdf', 'signal0', resol); 
print(typef, 'signal0', resol); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
plot(tt, x, 'k');

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
    
    % Erase x-axis :P
    plot([c, a]/fs, [0, 0], 'color', 'w', ...
            'linestyle', '-', ...
            'linewidth', 2.0);    
	
    % Draw first segment limit
    plot([a, a]/fs, [xmin, xmax], 'color', 'r', ...
            'linewidth', 3.0, ...
            'linestyle', '--');
    
    % Draw second segment limit
    plot([b, b]/fs, [xmin, xmax], 'color', 'r', ...
            'linewidth', 3.0, ...
            'linestyle', '--');
    
    c = b;    
    
    j = j + 1;
    zz(pp: (pp+pq)) = z;
    
    pp = pp + pq;
end
box off;
set(gca, 'FontSize', tit_fs);
ylabel('Amplitud', 'FontSize', lbl_fs);
xlabel('Tiempo (s)', 'FontSize', lbl_fs);

% Print figure to pdf and eps files    
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 3]);
set(gcf, 'PaperPosition', [0 0 8 3]);

print('-dpdf', 'signal1', resol); 
print(typef, 'signal1', resol); 

zz = zz(1:pp);

wavwrite(zz, fs, arch2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

min_x = 0;
max_x = length(zz)/fs;

max_y = max(ceil(abs(10*[min(zz), max(zz)])) / 10);
min_y = -max_y;

zz = [zz; zeros(length(x) - length(zz), 1)];
tt = (1:length(zz))/fs;
plot(tt, zz, 'k');

box off;
set(gca, 'FontSize', tit_fs);
ylabel('Amplitud', 'FontSize', lbl_fs);
xlabel('Tiempo (s)', 'FontSize', lbl_fs);

% Print figure to pdf and eps files    
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 3]);
set(gcf, 'PaperPosition', [0 0 8 3]);

print('-dpdf', 'signal2', resol); 
print(typef, 'signal2', resol); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc;

close all;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cc = melcepst(zz);
% cc = flipud(cc');
% subplot(3, 1, 2);
% imagesc(cc);


% [ka kb kc] = kmeans(cc', 300);
% hold on;
% subplot(3, 1, 3);

% imagesc([min_x, max_x], [min_y/4, max_y/4], kc');
% imagesc([min_x, max_x], [min_y/4, max_y/4], kc');
% plot(tt, zz, 'k');

% imagesc([min_x, max_x], [-0.01, 0.01], kc');
% plot([min_x, max_x], [min_y, max_y], tt, zz, 'k');

% hold off;
% colormap("summer");

