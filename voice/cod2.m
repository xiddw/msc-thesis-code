function [] = cod2(input, ppath)
    tic;

    ext = '.wav';

    if nargin < 2
        ppath = 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\voice\aud\';
    end
    
    wfile = strcat(ppath, input, ext)

    [x fs bps] = wavread(wfile);
%{
    xmin = min(x);
    xmax = max(x);

    dur = length(x)/fs;

    ws = 0.1 * fs;
    ws = ws + (mod(ws, 2) == 1);
    gx = ones(ws, 1)';

    xx = fftconv(abs(x), gx)';
    xx = xx(ws/2:end-(ws/2));
    tol_silence = quantile(xx, 0.25);
    pos_silence = find(xx < tol_silence);
    num_silence = length(pos_silence)-1;
    %{
    %}

    zz = x;
    %%{
    min_x = 0;
    max_x = length(zz)/fs;

    max_y = max(ceil(abs(10*[min(zz), max(zz)])) / 10);
    min_y = -max_y;
    %%}

    % tt = (1:length(zz))/fs;
      %}

    cc = melcepst(x, fs, 'W', 19, 30, 20, 10);
    cc = flipud(cc');
    csvwrite(strcat(ppath, input, '_ceps.csv'), cc');

    toc
end