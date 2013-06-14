function [] = cod2(input, ppath)
    tic;

    ext = '.wav';

    if nargin < 2
        ppath = 'E:\ESCUELA\CIMAT\4 Semestre\ST2\prog\voice\aud\';
    end
    
    wfile = strcat(ppath, input, ext)

    [x fs bps] = wavread(wfile);

    cc = melcepst(x, fs, 'W', 19, 30, 20, 10);
    cc = flipud(cc');
    csvwrite(strcat(ppath, input, '_ceps.csv'), cc');

    toc
end