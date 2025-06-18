% main_with_functions.m
clc;
clear all;

snrr = -2:1:30;
nb = 1593; % number of symbols
N = 10;
pilot = 1 + 1i; % pilot symbol

% 16-QAM constellation
inphase = [0.9 0.9 0.9 0.9 0.3 0.3 0.3 0.3];
quadr = [-0.9 -0.3 0.3 0.9 -0.9 -0.3 0.3 0.9];
inphase = [inphase; -inphase]; inphase = inphase(:);
quadr = [quadr; quadr]; quadr = quadr(:);
const = inphase + 1j * quadr;
M = length(const);

b = randi([0 M-1], nb, 1); % random symbols
tx = genqammod(b, const);

% Insert pilots
txp = addpilot(tx, nb, pilot, N);
len = length(txp);

% Rayleigh fading channel
ct = rayleighfading(len);
twn = txp .* ct; % apply channel

% Add noise
t = awgn(twn, 30, 'measured', 'db');

% Receiver: extract pilots
[rx, rxp] = extractpilot(t, N, nb);

% Channel estimation
csi = rxp / pilot;
chnstinf = interpft(csi, nb); % FFT interpolation

% Spline, linear, pchip interpolation
t1 = 1:1:1770;
t2 = 1:10:1770;
m = 1;
for i = 1:10:1770
    t3(m:m+8) = t1(i+1:i+9);
    m = m + 9;
end

chnstinf2 = interp1(t2, csi, t3, 'spline');
chnstinf3 = interp1(t2, csi, t3, 'linear');
chnstinf4 = interp1(t2, csi, t3, 'pchip');

% Zero Forcing Equalizer
for i = 1:nb
    RX(i) = rx(i) / chnstinf(i);
    RX2(i) = rx(i) / chnstinf2(i);
    RX3(i) = rx(i) / chnstinf3(i);
    RX4(i) = rx(i) / chnstinf4(i);
end

[ch, chp] = extpltcoef(ct, N, nb);

% BER simulation loop
for l = 1:length(snrr)
    tt = awgn(twn, snrr(l), 'measured', 'db');
    [rxt, rxpt] = extractpilot(tt, N, nb);
    csit = rxpt / pilot;

    chnstinft = interpft(csit, nb);
    t1t = 1:1:1770;
    t2t = 1:10:1770;
    mt = 1;
    for i = 1:10:1770
        t3t(mt:mt+8) = t1t(i+1:i+9);
        mt = mt + 9;
    end
    chnstinf2t = interp1(t2t, csit, t3t, 'spline');
    chnstinf3t = interp1(t2t, csit, t3t, 'linear');
    chnstinf4t = interp1(t2t, csit, t3t, 'pchip');

    for i = 1:nb
        RXt(i) = rxt(i) / chnstinft(i);
        RX2t(i) = rxt(i) / chnstinf2t(i);
        RX3t(i) = rxt(i) / chnstinf3t(i);
        RX4t(i) = rxt(i) / chnstinf4t(i);
    end

    rt = genqamdemod(RXt, const);
    r2t = genqamdemod(RX2t, const);
    r3t = genqamdemod(RX3t, const);
    r4t = genqamdemod(RX4t, const);

    [~, rate1(l)] = biterr(b.', rt);
    [~, rate2(l)] = biterr(b.', r2t);
    [~, rate3(l)] = biterr(b.', r3t);
    [~, rate4(l)] = biterr(b.', r4t);
end

% BER plot
figure;
semilogy(snrr, rate1, 'b-', snrr, rate2, 'r-', snrr, rate3, 'k-', snrr, rate4, 'g-');
legend('fft', 'cubic spline', 'linear', 'cubic');
title('BER curves for different interpolation techniques');
xlabel('SNR in dB');
ylabel('BER');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCIONES AUXILIARES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [txp] = addpilot(tx, nb, pilot, N)
    m = 1;
    txp = [];
    for i = 1:N-1:nb
        if i + N - 2 <= nb
            txp(m) = pilot;
            txp(m+1 : m+N-1) = tx(i : i+N-2);
            m = m + N;
        else
            txp(m) = pilot;
            txp(m+1 : m + nb - i) = tx(i : end);
            break;
        end
    end
end

function [rx, rxp] = extractpilot(t, N, nb)
    m = 1;
    rx = [];
    for i = 1:N:length(t)
        rxp(m) = t(i); % piloto
        rx = [rx t(i+1:i+N-1)]; % datos
        m = m + 1;
        if length(rx) >= nb
            break;
        end
    end
    rx = rx(1:nb); % ajustar longitud
end

function [ch, chp] = extpltcoef(ct, N, nb)
    m = 1;
    ch = [];
    for i = 1:N:length(ct)
        chp(m) = ct(i); % piloto
        ch = [ch ct(i+1:i+N-1)]; % datos
        m = m + 1;
        if length(ch) >= nb
            break;
        end
    end
    ch = ch(1:nb); % ajustar longitud
end

function [c] = rayleighfading(m)
    N = 32; % número de trayectorias
    fmax = 100; % frecuencia Doppler máxima
    A = 1; % amplitud
    f = 10000; % frecuencia de muestreo
    t = 0:1/f:(m/f - 1/f); % tiempo de muestreo
    ct = zeros(1, m);
    ph = 2 * pi * rand(1, N);
    theta = 2 * pi * rand(1, N);
    fd = fmax * cos(theta); % desplazamientos Doppler
    for n = 1:m
        for i = 1:N
            ct(n) = ct(n) + A * exp(1j * (2 * pi * fd(i) * t(n) + ph(i)));
        end
    end
    c = ct / sqrt(N); % normalización
end
