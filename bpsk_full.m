clc
clear all
close all

% Parámetros del sistema
nFFT = 64;               % Tamaño FFT
nDSC = 52;               % Número de subportadoras de datos
nBitPerSym = 52;         % Bits por símbolo OFDM (BPSK)
nSym = 1e4;              % Número de símbolos OFDM
EbN0dB = 0:2:30;         % Rango de Eb/N0
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % Conversión a Es/N0

% Parámetros del canal móvil
fc = 2.4e9;              % Frecuencia portadora
lambda = 3e8/fc;         % Longitud de onda
v = (30)/3.6;            % Velocidad del móvil (50 km/h -> m/s)
fmax = v/lambda;         % Frecuencia Doppler máxima
nTap = 10;               % Número de taps del canal
N = 50;                  % Número de paths para el modelo de movilidad

% Prealocación para almacenar BER
simBer = zeros(1, length(EbN0dB));
theoryBerAWGN = 0.5*erfc(sqrt(10.^(EbN0dB/10)));  % <-- Aquí estaba el error
theoryBerRayleigh = 0.5*(1-sqrt(10.^(EbN0dB/10)./(10.^(EbN0dB/10)+1)));

for ii = 1:length(EbN0dB)
    % Transmisor
    ipBit = rand(1, nBitPerSym*nSym) > 0.5;  % Bits aleatorios
    ipMod = 2*ipBit-1;                       % Modulación BPSK
    ipMod = reshape(ipMod, nBitPerSym, nSym).';
    
    % Mapeo a subportadoras
    xF = [zeros(nSym,6) ipMod(:,1:nBitPerSym/2) zeros(nSym,1) ...
          ipMod(:,nBitPerSym/2+1:nBitPerSym) zeros(nSym,5)];
    
    % IFFT
    xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';
    
    % Prefijo cíclico
    xt = [xt(:,49:64) xt];
    
    % Modelado del canal con los 3 fenómenos
    % 1. Canal multitrayecto con movilidad (Jakes model)
    ht = zeros(nSym, nTap);
    for symIdx = 1:nSym
        % Generación del canal para cada símbolo con efecto Doppler
        H = zeros(1, nTap);
        for n = 1:N
            an = sqrt(1/N);
            thetan = 2*pi*rand;
            fDn = fmax * cos(2*pi*rand);
            t = symIdx/nSym;  % Tiempo normalizado
            phase = thetan - 2*pi*fDn*t;
            H = H + an*exp(1j*phase);
        end
        ht(symIdx,:) = H/sqrt(N);  % Normalización de potencia
    end
    
    % 2. Convolución con el canal (efecto multitrayecto)
    xt_conv = zeros(nSym, 80+nTap-1);
    for jj = 1:nSym
        xt_conv(jj,:) = conv(ht(jj,:), xt(jj,:));
    end
    
    % 3. Conversión a vector largo y adición de AWGN
    xt_conv = reshape(xt_conv.', 1, nSym*(80+nTap-1));
    nt = 1/sqrt(2)*[randn(1, nSym*(80+nTap-1)) + 1j*randn(1, nSym*(80+nTap-1))];
    yt = sqrt(80/64)*xt_conv + 10^(-EsN0dB(ii)/20)*nt;
    
    % Receptor
    yt = reshape(yt.', 80+nTap-1, nSym).';
    yt = yt(:,17:80);  % Eliminación del prefijo cíclico
    
    % FFT
    yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).';
    
    % Estimación y ecualización del canal (se asume estimación perfecta)
    hF = fftshift(fft(ht, 64, 2));
    yF = yF./hF;
    
    % Extracción de subportadoras de datos
    yMod = yF(:,[6+(1:nBitPerSym/2) 7+(nBitPerSym/2+1:nBitPerSym)]);
    
    % Demodulación BPSK
    ipModHat = sign(real(yMod));
    ipBitHat = (ipModHat+1)/2;
    ipBitHat = reshape(ipBitHat.', nBitPerSym*nSym, 1).';
    
    % Cálculo de errores
    nErr = sum(ipBitHat ~= ipBit);
    simBer(ii) = nErr/(nBitPerSym*nSym);
end

% Gráficos
figure;
semilogy(EbN0dB, theoryBerAWGN, 'b-', 'LineWidth', 2);
hold on;
semilogy(EbN0dB, theoryBerRayleigh, 'r--', 'LineWidth', 2);
semilogy(EbN0dB, simBer, 'ks-', 'LineWidth', 2);
grid on;
axis([0 30 1e-5 1]);
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
legend('AWGN Teórico', 'Rayleigh Teórico', 'Simulación (AWGN+Rayleigh+Movilidad)', 'Location', 'southwest');
title('BER para OFDM-BPSK con AWGN, Rayleigh y Movilidad');

%% Visualización de la respuesta del canal
figure;
plot(abs(ht(1,:))), title('Respuesta al impulso del canal (1er símbolo)');
xlabel('Tap'); ylabel('Magnitud');

figure;
plot(20*log10(abs(hF(1,:)))), title('Respuesta en frecuencia del canal (1er símbolo)');
xlabel('Subportadora'); ylabel('Magnitud (dB)');