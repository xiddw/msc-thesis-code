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

font = 12;
tit_fs = 14;

% EXAMPLE Simple demo of the MFCC function usage.
%
%   This script is a step by step walk-through of computation of the
%   mel frequency cepstral coefficients (MFCCs) from a speech signal
%   using the MFCC routine.
%
%   See also MFCC, COMPARE.

%   Author: Kamil Wojcicki, September 2011


    % Clean-up MATLAB's environment
    %  clear all; close all; clc;  

    
    % Define variables
    Tw = 200;                % analysis frame duration (ms)
    Ts = 100;                % analysis frame shift (ms)
    alpha = 0.97;           % preemphasis coefficient
    M = 30;                 % number of filterbank channels 
    C = 19;                 % number of cepstral coefficients
    L = 22;                 % cepstral sine lifter parameter
    LF = 300;               % lower frequency limit (Hz)
    HF = 3700;              % upper frequency limit (Hz)
    wav_file  = 'voice\aud\amorosos_fin.wav';  % input audio filename    
    mfcc_file = 'voice\aud\amorosos_fin.csv';
    
    disp(strcat('Input: ', wav_file));
    disp(strcat('Output: ', mfcc_file));

    % Read speech samples, sampling rate and precision from file
    [ speech, fs, nbits ] = wavread( wav_file );
    % speech = speech(1:int64(length(speech )/10));

    % Feature extraction (feature vectors as columns)
    [ MFCCs, FBEs, frames, H, F ] = ...
                    mfcc( speech, fs, Tw, Ts, alpha, @hamming, [LF HF], M, C+1, L );
           
	figure;
        
    set(0, 'DefaultAxesFontSize', font)

    HH = H';
	plot(F(1:1000), HH(1:1000, :), 'LineWidth', 2);
    xlabel('Frecuencia (Hz)', 'FontSize', tit_fs);
    ylabel('Magnitud', 'FontSize', tit_fs);
    
    box off;
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 3]);
    set(gcf, 'PaperPosition', [0 0 8 3]);

    arch = 'mfcc_bankfilter';
    print(typef, arch, resol);     
    print('-dpdf', arch, resol);
      
    figure;    
    csvwrite(mfcc_file, MFCCs(2:end, :));


    % Generate data needed for plotting 
    [ Nw, NF ] = size( frames );                % frame length and number of frames
    time_frames = (0:(NF-1))*Ts*0.001+0.5*Nw/fs;  % time vector (s) for frames 
    time = (0:(length(speech)-1))/fs;           % time vector (s) for signal samples 
    logFBEs = 20*log10( FBEs );                 % compute log FBEs for plotting
    logFBEs_floor = max(logFBEs(:))-50;         % get logFBE floor 50 dB below max
    logFBEs( logFBEs<logFBEs_floor ) = logFBEs_floor; % limit logFBE dynamic range

    sp = [];
    sl = [];
    st = [];
    
    num = 400;
    idx = 1;

    % Generate plots
    %{
    figure('Position', [30 30 800 600], 'PaperPositionMode', 'auto', ... 
              'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' ); 
    %}
    
    mfilename = 'mfcc_result';

    % sp(1) = subplot( num + 10 + idx );
    % idx = idx + 1;
    plot(time, speech, 'k');
    box off;
    ylim([min(speech)-0.5 max(speech)+0.5]);
    xlim([min(time_frames) max(time_frames)]);
    xlabel('Tiempo (s)', 'FontSize', tit_fs); 
    ylabel('Amplitud', 'FontSize', tit_fs); 
    % st(1) = title('Nueva señal de audio'); 
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 3]);
    set(gcf, 'PaperPosition', [0 0 8 3]);
    
    print('-dpdf', sprintf('%s1.pdf', mfilename), resol); 
    print(typef, sprintf('%s1.eps', mfilename), resol);     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure; 

    % sp(2) = subplot( num + 10 + idx );
    % idx = idx + 1;
    imagesc( time_frames, (1:C), logFBEs); 
    colorbar();
    axis('xy');
    xlim([min(time_frames) max(time_frames)]);
    xlabel('Time (s)', 'FontSize', tit_fs); 
    ylabel('Num. de Coef.', 'FontSize', tit_fs); 
    % title( 'Respuesta al banco de filtros Mel (log)', 'FontSize', 16); 
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 3]);
    set(gcf, 'PaperPosition', [0 0 8 3]);
    
    print('-dpdf', sprintf('%s2.pdf', mfilename), resol); 
    print(typef, sprintf('%s2.eps', mfilename), resol);     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;     
    % sp(3) = subplot( num + 10 + idx );
    % idx = idx + 1;
    imagesc(time_frames, (1:C), MFCCs(2:end,:)); % HTK's TARGETKIND: MFCC
    colorbar();
    axis('xy');
    xlim([ min(time_frames) max(time_frames) ]);
    xlabel('Tiempo (s)', 'FontSize', tit_fs); 
    ylabel('Num. de Coef.', 'FontSize', tit_fs);
    % st(3) = title( 'Mel frequency cepstrum coefficient' );
    % Print figure to pdf and eps files    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 3]);
    set(gcf, 'PaperPosition', [0 0 8 3]);
    
    print('-dpdf', sprintf('%s3.pdf', mfilename), resol); 
    print(typef, sprintf('%s3.eps', mfilename), resol);     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure; 
    % sp(4) = subplot( num + 10 + idx );
    % idx = idx + 1;
    aa = kmeans(MFCCs(2:end, :)', 25);
    imagesc(time_frames, (1:C), aa' );
    colorbar();
    % st(4) = title( 'Coeficientes agrupados con k-means' );
    ylabel('Clusters', 'FontSize', tit_fs);  
    xlabel('Tiempo (s)', 'FontSize', tit_fs);    
    box off;
    set(gca, 'yticklabel', []);

    % Print figure to pdf and eps files    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 3]);
    set(gcf, 'PaperPosition', [0 0 8 3]);
    
    print('-dpdf', sprintf('%s4.pdf', mfilename), resol); 
    print(typef, sprintf('%s4.eps', mfilename), resol); 
    
    close all;
    
% EOF