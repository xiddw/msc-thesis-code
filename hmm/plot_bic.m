%%%%%%%%%%%%%$$$$$$$$$$$$$$
%%%%%%%%%%%%%$$$$$$$$$$$$$$

resol = '-r100';
typef = '-depsc';

set(0, 'DefaultAxesFontSize', 14);
tit_fs = 18;

if ~exist('K', 'var') && exist('orig', 'var')
    K = max(orig.obs);
end

if ~exist('N', 'var') && exist('orig', 'var')
    N = max(orig.hid);
end 

if ~exist('T', 'var') && exist('orig', 'var')
    T = length(orig.obs);
end

lll = (-10e2:1e2:60e2);
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
figure; 
% set(0, 'DefaultAxesFontSize', 28);

h = surfc(lll, seq_offs+(1:hh), bb, 'EdgeColor', 'white'); 
sx = lll(:);
sy = seq_offs+(1:hh);
sz = bb(:);
stez = range(sz)/5;
minz = str2double(sprintf('%.1f', min(sz)));
maxz = str2double(sprintf('%.1f', max(sz)));
stez = str2double(sprintf('%.1f', stez));

set(gca, 'XTick', min(sx):range(sx)/5:max(sx));
set(gca, 'YTick', 0:2:10);
set(gca, 'ZTick', minz:stez:maxz);
%rotate(h, [1 1 0], -0.5);
colormap(ametrine);

%shading interp;
xlabel('A', 'FontSize', tit_fs);
ylabel('B', 'FontSize', tit_fs);
zlabel('C', 'FontSize', tit_fs);
    
%xlabel('lambda', 'FontSize', tit_fs);
%ylabel('Modelo seleccionado', 'FontSize', tit_fs);
%zlabel('log-verosimilitud', 'FontSize', tit_fs);

%title(sprintf('Superficie de curvas BIC para distintos valores de lambda'));
if exist('archivo', 'var') 
    name = archivo;
    if exist('farchivo', 'var') 
        name = farchivo;
    end
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 5]);
    set(gcf, 'PaperPosition', [0 0 8 5]);

    arch1 = strcat(name, 'bic1');
    print(typef, arch1, resol);
end


%%% TEMP %%%
%%% xt = str2num(get(gca, 'XTickLabel'));

x = 1:length(lll);
xx = lll;
y = seq_offs+(1:hh);
z = bb;

st = max(diff(xx));
yy = st*y;

[px, py] = gradient(z, 2, 2);

figure;
pp = (abs(px) + abs(py)); % sqrt(abs(px.* py)); % min(abs(px), abs(py));
tot = sum(pp, 1);
[~, mp] = min(tot);
ar = cool(3);
%%{
%ar = [1, 0.8, 1; 1, 1, 1; 0.9, 1, 1];
n = 256; % size of new color map
m = size(ar,1);
t0 = linspace(0,1,m)';
t = linspace(0,1,n)';
r = interp1(t0,ar(:,1),t);
g = interp1(t0,ar(:,2),t);
b = interp1(t0,ar(:,3),t);
rar = [r,g,b];
%}

%imagesc(xx, yy, pp);
%colormap(ametrine(100));
%colorbar();

set(gca,'YDir','normal')
box off;
freezeColors 
hold on;

pc = contour(xx, yy, z, 30, 'LineWidth', 2); 
colormap(ametrine(100, 'invert', 1));
caxis auto

xlabel('lambda', 'FontSize', tit_fs)
ylabel('Modelo seleccionado', 'FontSize', tit_fs)

%figure;
quiver(xx, yy, px, py, 'k', 'LineWidth', 2);

yt = str2num(get(gca, 'YTickLabel'));
set(gca,'YTickLabel', yt/st);

plot([xx(mp), xx(mp)], [min(yy)-2, max(yy)+2], '--b', 'LineWidth', 3);

%tt = title(sprintf('Curva de nivel de superficie BIC para distintos valores de lambda'));
if exist('archivo', 'var') 
    name = archivo;
    if exist('farchivo', 'var') 
        name = farchivo;
    end
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [16 5]);
    set(gcf, 'PaperPosition', [0 0 16 5])
    
    arch2 = strcat(name, 'bic2');
    print(typef, arch2, resol);
end

%%
figure;
%set(0, 'DefaultAxesFontSize', 28);
mp = mp;
lambda = xx(mp);

ll = bb(:, mp);
[~, ml] = max(ll);

sp(1) = plot(y, ll, 'b');
xlabel('Modelo seleccionado', 'FontSize', tit_fs)
ylabel('log-verosimilitud', 'FontSize', tit_fs)

hold on;
sp(2) = plot(y, ll, 'og', 'MarkerFaceColor', 'g');
sp(3) = plot(y(ml), ll(ml), 'or', 'MarkerFaceColor', 'r');
% title('Selección de modelo con BIC');

set(gca, 'box', 'off');
set(sp(2:3),'MarkerSize', 7);
set(sp(3),'MarkerSize', 10);
set(sp, 'linewidth', 2)

l = legend(strcat('lambda= ', int2str(lambda)), 'Location', 'SouthEast');
set(l, 'Box', 'off');

if exist('archivo', 'var')
    name = archivo;
    if exist('farchivo', 'var') 
        name = farchivo;
    end
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 5]);
    set(gcf, 'PaperPosition', [0 0 8 5]);
    
    arch3 = strcat(name, 'bic3');
    print(typef, arch3, resol);
end

close all;
%{
print('-dpng', '', resol);
%}