clc
clear all
close all

% Parámetros del sistema
nFFT = 64;               % Tamaño FFT
nDSC = 52;               % Subportadoras de datos
nBitPerSym = 52;         % Bits por símbolo OFDM (BPSK)
nSym = 1e4;              % Número de símbolos OFDM
SNRdB = -2:1:30;         % Rango de SNRdB solicitado
EbN0dB = SNRdB;
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); 

% Canal móvil
fc = 2.4e9;              
lambda = 3e8/fc;         
v = (1)/3.6;             % Velocidad móvil (m/s)
fmax = v/lambda;         
nTap = 10;               
N = 1;                   % Número de paths

% Prealocación
simBer = zeros(1, length(SNRdB));
theoryBerAWGN = 0.5*erfc(sqrt(10.^(EbN0dB/10)));
theoryBerRayleigh = 0.5*(1 - sqrt(10.^(EbN0dB/10) ./ (10.^(EbN0dB/10) + 1)));

for ii = 1:length(SNRdB)
    % Transmisor
    ipBit = rand(1, nBitPerSym*nSym) > 0.5;
    ipMod = 2*ipBit - 1;                         % BPSK
    ipMod = reshape(ipMod, nBitPerSym, nSym).';

    % Mapeo a subportadoras
    xF = [zeros(nSym,6) ipMod(:,1:26) zeros(nSym,1) ipMod(:,27:52) zeros(nSym,5)];
    
    % IFFT
    xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';
    
    % Prefijo cíclico
    xt = [xt(:,49:64) xt];

    %% --- Visualización de la constelación transmitida ---
    if ii == 1
        figure;
        plot(real(ipMod(:)), imag(ipMod(:)), 'b.');
        title('Constelación Transmitida (BPSK)');
        xlabel('Re'); ylabel('Im'); grid on; axis equal;
    end

    %% --- AWGN puro: constelación sin canal multipath ---
    if ismember(SNRdB(ii), [-5 0 10 30])
        xt_awgn = reshape(xt.', 1, []);
        nt_awgn = 1/sqrt(2)*(randn(1, length(xt_awgn)) + 1j*randn(1, length(xt_awgn)));
        yt_awgn = sqrt(80/64)*xt_awgn + 10^(-EsN0dB(ii)/20)*nt_awgn;

        yt_awgn = reshape(yt_awgn.', 80, nSym).';
        yt_awgn = yt_awgn(:,17:80);
        yF_awgn = (sqrt(nDSC)/nFFT)*fftshift(fft(yt_awgn.')).';
        yMod_awgn = yF_awgn(:, [7:32, 34:59]);

        figure;
        idx = randperm(numel(yMod_awgn), 1000);
        plot(real(yMod_awgn(idx)), imag(yMod_awgn(idx)), '.');
        title(['Constelación con AWGN (SNR = ' num2str(SNRdB(ii)) ' dB)']);
        xlabel('Re'); ylabel('Im'); grid on; axis equal;
    end

    %% --- Canal móvil con movilidad y multitrayecto ---
    ht = zeros(nSym, nTap);
    for symIdx = 1:nSym
        H = zeros(1, nTap);
        for n = 1:N
            an = sqrt(1/N);
            thetan = 2*pi*rand;
            fDn = fmax * cos(2*pi*rand);
            t = symIdx / nSym;
            phase = thetan - 2*pi*fDn*t;
            H = H + an * exp(1j*phase);
        end
        ht(symIdx,:) = H / sqrt(N);
    end

    % Convolución con canal
    xt_conv = zeros(nSym, 80+nTap-1);
    for jj = 1:nSym
        xt_conv(jj,:) = conv(ht(jj,:), xt(jj,:));
    end

    % Canal + AWGN
    xt_conv = reshape(xt_conv.', 1, []);
    nt = 1/sqrt(2)*(randn(1, length(xt_conv)) + 1j*randn(1, length(xt_conv)));
    yt = sqrt(80/64)*xt_conv + 10^(-EsN0dB(ii)/20)*nt;

    % Receptor
    yt = reshape(yt.', 80+nTap-1, nSym).';
    yt = yt(:,17:80);  % Quitar CP

    yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).';

    %% --- Constelación antes de la ecualización ---
    if ismember(SNRdB(ii), [-5 0 10 30])
        yMod_preEq = yF(:, [7:32, 34:59]);
        figure;
        idx = randperm(numel(yMod_preEq), 1000);
        plot(real(yMod_preEq(idx)), imag(yMod_preEq(idx)), '.');
        title(['Canal+AWGN sin ecualizar (SNR = ' num2str(SNRdB(ii)) ' dB)']);
        xlabel('Re'); ylabel('Im'); grid on; axis equal;
    end

    %% --- Ecualización con CSI perfecto ---
    hF = fftshift(fft(ht, 64, 2));
    yF = yF ./ hF;

    yMod = yF(:, [7:32, 34:59]);

    %% --- Constelación después de ecualización ---
    if ismember(SNRdB(ii), [-5 0 10 30])
        figure;
        idx = randperm(numel(yMod), 1000);
        plot(real(yMod(idx)), imag(yMod(idx)), '.');
        title(['Canal+AWGN ecualizado (SNR = ' num2str(SNRdB(ii)) ' dB)']);
        xlabel('Re'); ylabel('Im'); grid on; axis equal;
    end

    % Demodulación y BER
    ipModHat = sign(real(yMod));
    ipBitHat = (ipModHat + 1)/2;
    ipBitHat = reshape(ipBitHat.', nBitPerSym*nSym, 1).';
    nErr = sum(ipBitHat ~= ipBit);
    simBer(ii) = nErr / (nBitPerSym*nSym);
end

%% --- Gráfica de BER ---
figure;
semilogy(SNRdB, theoryBerAWGN, 'b-', 'LineWidth', 2); hold on;
semilogy(SNRdB, theoryBerRayleigh, 'r--', 'LineWidth', 2);
semilogy(SNRdB, simBer, 'ks-', 'LineWidth', 2);
grid on; axis([-2 30 1e-5 1]);
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
legend('AWGN Teórico', 'Rayleigh Teórico', 'Simulación', 'Location', 'southwest');
title('BER para OFDM-BPSK con canal móvil');

%% --- Visualización del canal ---
figure;
plot(abs(ht(1,:))), title('Respuesta al impulso del canal');
xlabel('Tap'); ylabel('Magnitud');

figure;
plot(20*log10(abs(hF(1,:)))), title('Respuesta en frecuencia del canal');
xlabel('Subportadora'); ylabel('Magnitud (dB)');
