clc; 
clear all;
close all;
%% Task1
%%sampling
m=@(t)((5.*(cos(20*pi*t)))+ (4.*(sin(10*pi*t))));%Message signal
Fs = 100; %Sampling Frequency
t = 0:0.0001:0.1; % Time Grid
ts = 0:1/Fs:0.1; % Sampling Time Grid
ms = m(ts); % Sampled Signal
figure(1)
subplot(2,1,1), plot(t, m(t), 'linewidth', 1);
title('Sampling Process');
ylabel('Message signal, m(t)');
subplot(2,1,2), stem(ts, ms, '.', 'linewidth', 1.5);
ylabel('Sampled Message Signal, m_s(n)');
xlabel('Time (sec)')

%% Quantization (Rounding)
L = 256; b = log2(L); %QuantizationLevels and bits
delta = (max(ms)-min(ms))/(L-1); % Quantization Step
xq = round(ms/delta)*delta; % Quantized Signal
SQNR = 10*log10(mean(ms.^2)/mean((ms-xq).^2));%Signal-to-quantization-noise ratio
disp(['SQNR = ', num2str(SQNR), ' dB']);
figure(2)
subplot(2,1,1), stem(ts, ms, '.', 'linewidth', 1.5);
title('Quantization Process');
ylabel('Sampled Signal, x_s(n)');
subplot(2,1,2), stem(ts, xq, '.', 'linewidth', 1.5);
ylabel('Quantized Signal, x_q(n)');

%% PCM Encoding
xb = de2bi(round((xq/delta)+41)); %Assigning bits foreach Q-level may addl/2 
xbs = reshape(xb, [1, length(xq)*b]); %PCM Bitstream
figure(3)
stairs(xbs, 'linewidth', 1);
title('PCM Encoding');

ylabel('PCM Bitstream');
xlim([0, length(xbs)-1]);
ylim([-0.05, 1.25]);



%% Reconstruction
xbr = reshape(xbs, [length(xbs)/b, b]);
xqr = (bi2de(xbr) - L/2); % reconstructed quantized signal
xr = interp1(ts, xqr, t); % reconstructed message
figure(4)
subplot(2,1,1), stem(ts, xqr, '.', 'linewidth', 1.5);
title('Reconstructed Signal');
ylabel('Reconstructed x_q(n)');
subplot(2,1,2), plot(t, xr, 'linewidth', 1);
ylabel('Reconstructed Message Signal, x_r(t)');
xlabel('Time (sec)');
%% Task2 :: Modulation of the PCM bitstream using M-ary QAM modulation scheme by taking M=16
M=16;
k=log2(M);
dataInMatrix = reshape(xbs,length(xbs)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInMatrix);     % Convert to integers

figure(5); % Create new figure window.
stem(dataSymbolsIn(1:16));
title('Random Symbols');
xlabel('Symbol Index');
ylabel('Integer Value');


dataMod = qammod(dataSymbolsIn,M);  %bitsteam modulated to QAM
figure(6)
stem(dataMod,'linewidth',1);
title('Quadrature Amplitude Modulated signal')

%% Task 3:: Addition of additive white Gaussian noise of 10dB with the QAM modulated signal
numsamplespersymbol=1;

EbNo=10;%Additive white Guassian Noise
snr=EbNo+10*log10(numsamplespersymbol);%Signal to Noise Ratio

noisysignal=awgn(dataMod,snr,'measured');
sPlotFig = scatterplot(noisysignal,1,0,'g.');
hold on
figure(7)
scatterplot(dataMod,1,0,'k*',sPlotFig);
figure(8)
stem(noisysignal);
title('QAM signal after adding AWGN ')
%%TASK 4 (Demodulating the modulated qam signal)
M=M;
dataout=qamdemod(noisysignal,M);
figure(9)
stem(dataout);
title('demodulated QAM signal')

%% TASK 4:: Part 1: Recovering the PCM bitstream
dataoutMatrix= de2bi(dataout,log2(M));
dataoutbitstream=dataoutMatrix(:);
figure(10)
stairs(dataoutbitstream)
 title('Recovered PCM Bitstream')
%% TASK 4:: Part 2: Comparing with actual bitstreams and finding BER)
[numErrors,ber]=biterr(dataInMatrix,dataoutMatrix);
fprintf('\nThe bit error rate = %5.2e, based on %d errors\n', ...
    ber,numErrors)

%% Task5:: Recovery of the signal from recovered bitstream

xbrr = reshape(dataoutbitstream, [length(dataoutbitstream)/b, b]);
xqrr = (bi2de(xbrr) - L/2); % reconstructed quantized signal
xrr = interp1(ts, xqrr, t); % reconstructed message
figure(11)
subplot(2,1,1), stem(ts, xqrr, '.', 'linewidth', 1.5);
title('Reconstructed Signal');
ylabel('Reconstructed x_q(n)');
subplot(2,1,2), plot(t, xrr, 'linewidth', 1);
ylabel('Reconstructed Message Signal, x_r(t)');
xlabel('Time (sec)');
%%quantization
L = 256; b = log2(L); %QuantizationLevels and bits
delta = (max(ms)-min(ms))/(L-1); % Quantization Step
xqrr = round(ms/delta)*delta; % Quantized Signal
SQNR = 10*log10(mean(ms.^2)/mean((ms-xqrr).^2));
disp(['SQNR = ', num2str(SQNR), ' dB']);
figure(12)
subplot(2,1,1), stem(ts, ms, '.', 'linewidth', 1.5);
title('Quantization Process after recieving');
ylabel('Sampled Signal, x_s(n)');
subplot(2,1,2), stem(ts, xqrr, '.', 'linewidth', 1.5);
ylabel('Quantized Signal after recieving, x_q(n)');
%%%%done by Ayan Chowdhury
