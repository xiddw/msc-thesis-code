%%%%%%%%%%%%%$$$$$$$$$$$$$$
%%%%%%%%%%%%%$$$$$$$$$$$$$$

resol = '-r400';

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
set(0, 'DefaultAxesFontSize', 13)

surfc(lll, seq_offs+(1:hh), bb); colormap('cool');
xlabel('lambda')
ylabel('Modelo seleccionado')
zlabel('log-verosimilitud')
tt = title(sprintf('Superficie de curvas BIC para distintos valores de lambda'));
set(tt, 'FontSize', 16)

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
% ar = colormap('cool');
ar = [1, 0.8, 1; 1, 1, 1; 0.9, 1, 1];
% ar = ar(end:-1:1, :);

n = 256; % size of new color map
m = size(ar,1);
t0 = linspace(0,1,m)';
t = linspace(0,1,n)';
r = interp1(t0,ar(:,1),t);
g = interp1(t0,ar(:,2),t);
b = interp1(t0,ar(:,3),t);
rar = [r,g,b];

imagesc(xx, yy, pp);
colormap(rar)
xlabel('lambda')
ylabel('Modelo seleccionado')

set(gca,'YDir','normal')
box off;
freezeColors 
hold on;

pc = contour(xx, yy, z, 30, 'LineWidth', 2); 
plot([xx(mp), xx(mp)], [min(yy)-2, max(yy)+2], '--p', 'LineWidth', 3);
colormap('cool');
caxis auto

%figure;
quiver(xx, yy, px, py, 'k', 'LineWidth', 1);

yt = str2num(get(gca, 'YTickLabel'));
set(gca,'YTickLabel', yt/st);

tt = title(sprintf('Curva de nivel de superficie BIC para distintos valores de lambda'));
set(tt, 'FontSize', 16)

%%
figure;
mp = mp;
lambda = xx(mp);

ll = bb(:, mp);
[~, ml] = max(ll);

sp(1) = plot(y, ll, 'b');
xlabel('Modelo seleccionado')
ylabel('log-verosimilitud')

hold on;
sp(2) = plot(y, ll, 'or', 'MarkerFaceColor', 'r');
sp(3) = plot(y(ml), ll(ml), 'og', 'MarkerFaceColor', 'g');
t1 = title('Selección de modelo con BIC');

set(t1, 'FontSize', 16)
set(gca, 'FontSize', 12);
set(gca, 'box', 'off');
set(sp(2:3),'MarkerSize', 7);
set(sp, 'linewidth', 2)

legend(strcat('lambda= ', int2str(lambda)))

if exist('archivo', 'var') 
    arch1 = strcat(archivo, '_1');
    arch2 = strcat(archivo, '_2');
    arch3 = strcat(archivo, '_3');
    
    print('-dpng', arch1, resol); close;
    print('-dpng', arch2, resol); close;
    print('-dpng', arch3, resol); %close;
end
%{
print('-dpng', '', resol);
%}