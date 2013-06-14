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
    wav_file  = 'voice\aud\cats1f.wav';  % input audio filename    
    mfcc_file = 'voice\aud\cats1f_ceps.csv';
    
    disp(strcat('Input: ', wav_file));
    disp(strcat('Output: ', mfcc_file));

    % Read speech samples, sampling rate and precision from file
    [ speech, fs, nbits ] = wavread( wav_file );


    % Feature extraction (feature vectors as columns)
    [ MFCCs, FBEs, frames ] = ...
                    mfcc( speech, fs, Tw, Ts, alpha, @hamming, [LF HF], M, C+1, L );
                
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

    sp(1) = subplot( num + 10 + idx );
    idx = idx + 1;
    plot( time, speech, 'k' );
    %xlim( [ min(time_frames) max(time_frames) ] );
    %xlabel( 'Tiempo (s)' ); 
    sl(1) = ylabel( 'Amplitud' ); 
    %st(1) = title( 'Nueva señal de audio'); 

    sp(2) = subplot( num + 10 + idx );
    idx = idx + 1;
    imagesc( time_frames, [1:C], logFBEs ); 
    axis( 'xy' );
    xlim( [ min(time_frames) max(time_frames) ] );
    %xlabel( 'Time (s)' ); 
    sl(2) = ylabel( 'Num. de Coef.' ); 
    %st(2) = title( 'Respuesta al banco de filtros Mel (log)'); 

    sp(3) = subplot( num + 10 + idx );
    idx = idx + 1;
    imagesc( time_frames, [1:C], MFCCs(2:end,:) ); % HTK's TARGETKIND: MFCC
    %imagesc( time_frames, [1:C+1], MFCCs );       % HTK's TARGETKIND: MFCC_0
    axis( 'xy' );
    xlim( [ min(time_frames) max(time_frames) ] );
    %xlabel( 'Time (s)' ); 
    sl(3) = ylabel( 'Num. de Coef.' );
    %st(3) = title( 'Mel frequency cepstrum coefficient' );
    
    sp(4) = ...
    subplot( num + 10 + idx );
    idx = idx + 1;
    aa = kmeans(MFCCs(2:end, :)', 50);
    imagesc(time_frames, [1:C], aa' );
    %st(4) = ...
    %title( 'Coeficientes agrupados con k-means' )
    sl(4) = ylabel( 'Clusters' );  
    sl(5) = xlabel('Tiempo (s)');
    
    set(sp, 'FontSize', 12)
    set(sp, 'box', 'off')
    set(sp, 'color', 'white')
    set(sp, 'linewidth', 1)
    
    %set(sp(1:3), 'xticklabel', []);
    
    %set(sp(4), 'xticklabel', []);
    set(sp(4), 'yticklabel', []);

    set(sl, 'linewidth', 1.1)
    %set(st, 'FontSize', 16)
    set(sl, 'FontSize', 12)
    % Set color map to grayscale
    %colormap( 1-colormap('jet') ); 

    % Print figure to pdf and png files
    print('-dpdf', sprintf('%s.pdf', mfilename)); 
    print('-dpng', sprintf('%s.png', mfilename)); 
    
    
% EOF