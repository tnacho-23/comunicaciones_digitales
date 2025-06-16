clc; clear; close all;

% ----------- Parámetros  ----------
nBits = 1e5;                      % bits a transmitir
bps = 2;                          % Bits por símbolo (QPSK)
nDSC = 52;                        % Subportadoras útiles
nFFT = 64;                        % Tamaño de la FFT
cpLen = 16;                       % Largo del prefijo cíclico
nSym = ceil(nBits / (nDSC * bps)); % Símbolos OFDM necesarios

% ----------- Canal multipath Rayleigh -----------------------
L = 5;                     % Número de trayectorias
Ts = 1e-6;                 % Tiempo de muestreo
f_c = 700e6;               % Portadora central
v = 30 / 3.6;              % Velocidad (m/s)
c = 3e8;                   % Velocidad de la luz
fd = (v * f_c) / c;        % Doppler máximo

% ----------- SNRs a simular ---------------------------------
EbN0dB = -2:1:30;
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(nFFT / (nFFT + cpLen));
BER = zeros(size(EbN0dB));

% ----------- Mapeo QPSK (Gray) ------------------------------
mod_map = [1+1j, -1+1j, -1-1j, 1-1j]/sqrt(2); % Gray coding

% ----------- Generación de datos ----------------------------
ipBits = randi([0 1], nSym * nDSC * bps, 1);      % bits aleatorios
ipSymbols = reshape(ipBits, bps, []).';
ipDec = bi2de(ipSymbols, 'left-msb') + 1;
modData = mod_map(ipDec).';

% --- Constelación QPSK transmitida ---
figure;
scatter(real(modData), imag(modData), '.');
axis square; grid on;
title('Constelación QPSK transmitida (antes del canal)');
xlabel('Re'); ylabel('Im');

% ----------- Agrupar datos en símbolos OFDM -----------------
modData = reshape(modData, nDSC, nSym).';
xF = zeros(nSym, nFFT);
xF(:,[7:32 34:59]) = modData;  % asignar datos

% ----------- IFFT + prefijo cíclico -------------------------
xt = ifft(ifftshift(xF, 2), nFFT, 2) * sqrt(nFFT);
xt_cp = [xt(:,end-cpLen+1:end), xt];
tx = reshape(xt_cp.', [], 1);

% ----------- Simulación para cada SNR -----------------------
for ii = 1:length(EsN0dB)
    % Canal Rayleigh
    h = (1/sqrt(2*L)) * (randn(nSym, L) + 1j*randn(nSym, L));
    rx = zeros(size(tx));
    symLen = nFFT + cpLen;

    for n = 1:nSym
        h_sym = h(n,:);
        x_block = xt_cp(n,:);
        r_block = conv(x_block, h_sym, 'same');
        rx((n-1)*symLen+1:n*symLen) = r_block;
    end

    % Ruido AWGN
    noise = 1/sqrt(2)*(randn(size(rx)) + 1j*randn(size(rx)));
    rx_noisy = rx + 10^(-EsN0dB(ii)/20)*noise;

    % Receptor OFDM
    rx_blocks = reshape(rx_noisy, symLen, []).';
    rx_noCP = rx_blocks(:, cpLen+1:end);
    Y = fftshift(fft(rx_noCP, nFFT, 2), 2) / sqrt(nFFT);

    % Constelación con ruido (antes de equalizar)
    Y_no_eq = Y(:, [7:32 34:59]);
    if ismember(EbN0dB(ii), [-2 0 10 30])
        figure;
        scatter(real(Y_no_eq(:)), imag(Y_no_eq(:)), '.');
        axis square; grid on;
        title(['Constelación con ruido (sin equalizar), SNR = ', num2str(EbN0dB(ii)), ' dB']);
        xlabel('Re'); ylabel('Im');
    end

    % Perfect CSI y equalización ZF
    Y_used = Y(:,[7:32 34:59]);
    H_used = fftshift(fft(h, nFFT, 2), 2);
    H_used = H_used(:,[7:32 34:59]);
    eqSym = Y_used ./ H_used;

    % Constelación después de ecualización (ya existía)
    if ismember(EbN0dB(ii), [-2 0 10 30])
        figure;
        scatter(real(eqSym(:)), imag(eqSym(:)), '.');
        axis square; grid on;
        title(['Constelación QPSK ecualizada, SNR = ', num2str(EbN0dB(ii)), ' dB']);
        xlabel('Re'); ylabel('Im');
    end

    % Demodulación QPSK
    demod_real = real(eqSym) > 0;
    demod_imag = imag(eqSym) > 0;
    ipBitsHat = [demod_real(:) demod_imag(:)].';
    ipBitsHat = ipBitsHat(:);          % columna N×1
    ipRef     = ipBits(1:length(ipBitsHat));

    % BER
    BER(ii) = sum(ipBitsHat ~= ipRef) / length(ipRef);
end

% ----------- Curva teórica y gráfica final ------------------
EbN0lin = 10.^(EbN0dB/10);
theoryBer = 0.5*(1 - sqrt(EbN0lin ./ (1 + EbN0lin)));

figure;
semilogy(EbN0dB, BER, 'r-o', 'LineWidth', 2); hold on;
semilogy(EbN0dB, theoryBer, 'b--', 'LineWidth', 2);
xlabel('Eb/N0 [dB]'); ylabel('BER');
title('BER para OFDM-QPSK en canal Rayleigh con ZF (Perfect CSI)');
legend('Simulación', 'Teórico Rayleigh');
grid on;
