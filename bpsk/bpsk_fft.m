clc; clear; close all;

%% -------- Parámetros del sistema ---------------------------------------
nFFT        = 64;          % Tamaño FFT
nDSC        = 52;          % Subportadoras de datos
nBitPerSym  = nDSC;        % Bits por símbolo OFDM (BPSK)
nSym        = 1e4;         % Nº total de símbolos OFDM (datos + pilotos)
Npil        = 5;           % <-- NUEVO: separación de símbolos piloto
pilotIdx    = 1:Npil:nSym; % <-- NUEVO: posiciones de los símbolos piloto
dataIdx     = setdiff(1:nSym,pilotIdx);

SNRdB  = -2:1:30;
EbN0dB = SNRdB;
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80);

%% -------- Canal móvil (igual que antes) --------------------------------
fc  = 2.4e9;  lambda = 3e8/fc;
v   = 10/3.6; fmax   = v/lambda;
nTap = 10;  N = 5;

%% -------- Pre‑alocación -------------------------------------------------
simBer          = zeros(1,length(SNRdB));
thBerAWGN       = 0.5*erfc(sqrt(10.^(EbN0dB/10)));
thBerRayleigh   = 0.5*(1-sqrt(10.^(EbN0dB/10)./(10.^(EbN0dB/10)+1)));

%% -------- Piloto en frecuencia (todo BPSK conocido) --------------------
pilotPat = [zeros(1,6)  (-1).^(1:26)  0  (-1).^(27:52)  zeros(1,5)]; % |X_pilot|=1

for ii = 1:length(SNRdB)

    %% --------------------- TRANSMISOR ----------------------------------
    % Bits sólo para los símbolos de datos
    ipBit      = rand(1, nBitPerSym*length(dataIdx)) > 0.5;
    ipMod      = 2*ipBit - 1;                         % BPSK
    ipMod      = reshape(ipMod, nBitPerSym, []).';

    % Mapeo a subportadoras (todos los símbolos)
    xF = zeros(nSym,nFFT);
    
    % -- a) Símbolos piloto  (NUEVO) ------------------------------------
    xF(pilotIdx,:) = repmat(pilotPat,length(pilotIdx),1);
    
    % -- b) Símbolos de datos (MOD) -------------------------------------
    xF(dataIdx,:)  = [zeros(size(ipMod,1),6) ...
                      ipMod(:,1:26) zeros(size(ipMod,1),1) ...
                      ipMod(:,27:52) zeros(size(ipMod,1),5)];

    % IFFT y CP
    xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';
    xt = [xt(:,49:64) xt];          % añadir CP

    %% --------------------- CANAL ---------------------------------------
    ht = zeros(nSym,nTap);
    for s = 1:nSym
        H = 0;
        for n = 1:N
            an = sqrt(1/N);
            thetan = 2*pi*rand;
            fDn = fmax*cos(2*pi*rand);
            t = s/nSym;
            H = H + an*exp(1j*(thetan-2*pi*fDn*t));
        end
        ht(s,:) = H/sqrt(N);
    end

    % Convolución
    xt_conv = zeros(nSym,80+nTap-1);
    for s=1:nSym,  xt_conv(s,:) = conv(ht(s,:),xt(s,:)); end

    % Canal + AWGN
    xt_vec = reshape(xt_conv.',1,[]);
    nt     = (randn(size(xt_vec))+1j*randn(size(xt_vec)))/sqrt(2);
    yt_vec = sqrt(80/64)*xt_vec + 10^(-EsN0dB(ii)/20)*nt;

    %% --------------------- RECEPTOR ------------------------------------
    yt = reshape(yt_vec.',80+nTap-1,nSym).';
    yt = yt(:,17:80);                           % quitar CP
    yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).';

    % -------- Estimación de canal (NUEVO) -------------------------------
    hEst      = zeros(size(yF));
    % 1) estimación directa en símbolos piloto
    hEst(pilotIdx,:) = yF(pilotIdx,:) ./ xF(pilotIdx,:);
    % 2) interpolación lineal para los símbolos de datos
    for k = 1:nFFT
        hEst(:,k) = interp1(pilotIdx, hEst(pilotIdx,k), 1:nSym, 'linear', 'extrap');
    end

    % Igualación con hEst en vez de CSI perfecto
    yF_eq = yF ./ hEst;

    % Extracción de datos y demodulación
    dataSubs = [7:32 34:59];
    yMod     = yF_eq(dataIdx,dataSubs);
    ipHat    = sign(real(yMod));
    ipHat    = (ipHat+1)/2;
    ipHat    = reshape(ipHat.',1,[]);
    
    % BER
    nErr         = sum(ipHat ~= ipBit);
    simBer(ii)   = nErr/length(ipBit);

    % --- (Opcional) deja tus gráficas de constelación/BER aquí -----------
end

%% -------- Gráfica de BER ------------------------------------------------
semilogy(SNRdB,thBerAWGN,'b-', ...
         SNRdB,thBerRayleigh,'r--', ...
         SNRdB,simBer,'ks-','LineWidth',2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
legend('AWGN teórico','Rayleigh teórico','Simulación');
title(['BER con estimación FFT y pilotos cada ' num2str(Npil) ' símbolos']);