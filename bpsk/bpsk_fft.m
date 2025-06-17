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
pilotSpacing = 5;        % Pilotos cada 5 símbolos

% Canal móvil
fc = 2.4e9;              
lambda = 3e8/fc;         
v = 1/3.6;               % Velocidad móvil (m/s)
fmax = v/lambda;         
nTap = 10;               
N = 1;                   

% Prealocación
simBer = zeros(1, length(SNRdB));
theoryBerAWGN = 0.5*erfc(sqrt(10.^(EbN0dB/10)));
theoryBerRayleigh = 0.5*(1 - sqrt(10.^(EbN0dB/10) ./ (10.^(EbN0dB/10) + 1)));

for ii = 1:length(SNRdB)
    % Bits y pilotos
    ipBit = rand(1, nBitPerSym*nSym) > 0.5;
    ipMod = 2*ipBit - 1;                         
    ipMod = reshape(ipMod, nBitPerSym, nSym).';
    
    pilotSymbols = 2*randi([0 1], nSym/pilotSpacing, nBitPerSym) - 1;
    
    ipMod_pilots = ipMod;
    for p = 1:floor(nSym/pilotSpacing)
        ipMod_pilots(pilotSpacing*p,:) = pilotSymbols(p,:);
    end

    % Mapeo a subportadoras
    xF = [zeros(nSym,6) ipMod_pilots(:,1:26) zeros(nSym,1) ipMod_pilots(:,27:52) zeros(nSym,5)];
    
    xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';
    xt = [xt(:,49:64) xt];  % Añadir CP

    % Canal móvil
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

    % Convolución + canal + ruido
    xt_conv = zeros(nSym, 80+nTap-1);
    for jj = 1:nSym
        xt_conv(jj,:) = conv(ht(jj,:), xt(jj,:));
    end

    xt_conv = reshape(xt_conv.', 1, []);
    nt = 1/sqrt(2)*(randn(1, length(xt_conv)) + 1j*randn(1, length(xt_conv)));
    yt = sqrt(80/64)*xt_conv + 10^(-EsN0dB(ii)/20)*nt;

    yt = reshape(yt.', 80+nTap-1, nSym).';
    yt = yt(:,17:80);  % Quitar CP

    yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).';

    % Estimación de canal por pilotos
    estH_all = zeros(nSym, nFFT);
    for symIdx = 1:nSym
        if mod(symIdx, pilotSpacing) == 0
            txPilotF = xF(symIdx,:);
            estH = yF(symIdx,:) ./ txPilotF;
            estH_all(symIdx,:) = estH;
        end
    end

    % Interpolación del canal estimado
    pilotIdx = pilotSpacing:pilotSpacing:nSym;
    for k = 1:nFFT
        knownH = estH_all(pilotIdx,k);
        interpH = interp1(pilotIdx, knownH, 1:nSym, 'linear', 'extrap');
        estH_all(:,k) = interpH;
    end

    % Ecualización
    yF_eq = yF ./ estH_all;
    yMod = yF_eq(:, [7:32, 34:59]);

    % Visualizaciones
    if ismember(SNRdB(ii), [-2 0 10 30])
        % Constelación transmitida
        if ii == 1
            figure;
            plot(real(ipMod(:)), imag(ipMod(:)), 'b.');
            title('Constelación Transmitida (BPSK)');
            xlabel('Re'); ylabel('Im'); grid on; axis equal;
        end

        % Pre-ecualización
        figure;
        yPre = yF(:, [7:32, 34:59]);
        idx = randperm(numel(yPre), 1000);
        plot(real(yPre(idx)), imag(yPre(idx)), '.');
        title(['Antes de ecualizar (SNR = ' num2str(SNRdB(ii)) ' dB)']);
        xlabel('Re'); ylabel('Im'); grid on; axis equal;

        % Post-ecualización
        figure;
        idx = randperm(numel(yMod), 1000);
        plot(real(yMod(idx)), imag(yMod(idx)), '.');
        title(['Después de ecualizar (estimado) (SNR = ' num2str(SNRdB(ii)) ' dB)']);
        xlabel('Re'); ylabel('Im'); grid on; axis equal;
    end

    % Demodulación y BER
    ipModHat = sign(real(yMod));
    ipBitHat = (ipModHat + 1)/2;
    ipBitHat = reshape(ipBitHat.', nBitPerSym*nSym, 1).';
    nErr = sum(ipBitHat ~= ipBit);
    simBer(ii) = nErr / (nBitPerSym*nSym);

    % Guardar CSI perfecto solo para comparación
    if ii == 1
        hF_perfect = fftshift(fft(ht, 64, 2));
    end
end

% BER
figure;
semilogy(SNRdB, theoryBerAWGN, 'b-', 'LineWidth', 2); hold on;
semilogy(SNRdB, theoryBerRayleigh, 'r--', 'LineWidth', 2);
semilogy(SNRdB, simBer, 'ks-', 'LineWidth', 2);
grid on; axis([-2 30 1e-5 1]);
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
legend('AWGN Teórico', 'Rayleigh Teórico', 'Simulación', 'Location', 'southwest');
title('BER para OFDM-BPSK con canal móvil y estimación por FFT');

% Visualización del canal
figure;
plot(abs(ht(1,:)));
title('Respuesta al impulso del canal');
xlabel('Tap'); ylabel('Magnitud');

figure;
plot(20*log10(abs(hF_perfect(1,:))), 'b'); hold on;
plot(20*log10(abs(estH_all(1,:))), 'r--');
title('Respuesta en frecuencia del canal');
xlabel('Subportadora'); ylabel('Magnitud (dB)');
legend('CSI perfecto', 'Estimación FFT');
