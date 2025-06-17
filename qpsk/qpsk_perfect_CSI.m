clc; clear all; close all;

% Parámetros
num_bits = 1e5;
k = 2; % bits por símbolo para QPSK
num_symbols = num_bits / k;
EbN0_dB = -2:1:30;
L = 5; % número de trayectorias
v_kmh = 120; fc = 700e6;
v = v_kmh / 3.6;
lambda = 3e8 / fc;
fd_max = v / lambda;
num_runs = 1;

% Mapeo QPSK Gray
bit_map = [0 0; 0 1; 1 1; 1 0]; % Codificación Gray
mapping = [1+1j; -1+1j; -1-1j; 1-1j] / sqrt(2); % Símbolos QPSK

ber_total = zeros(1, length(EbN0_dB));

for idx = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(idx)/10);
    for run = 1:num_runs
        % Transmisión
        bits = randi([0 1], 1, num_bits);
        bit_pairs = reshape(bits, 2, []).';
        [~, indices] = ismember(bit_pairs, bit_map, 'rows');
        symbols = mapping(indices).';

        % Canal Rayleigh plano con L trayectorias
        t = linspace(0, 1, num_symbols);
        an = ones(1,L)/sqrt(L);
        thetan = 2*pi*rand(1,L);
        fDn = fd_max * cos(2*pi*rand(1,L));
        H = zeros(1,num_symbols);
        for l = 1:L
            H = H + an(l)*exp(1j*(thetan(l) - 2*pi*fDn(l)*t));
        end
        H(abs(H) < 1e-3) = 1e-3; % Evitar divisiones por 0

        % Canal con ruido AWGN usando la función awgn (SNR por símbolo)
        snr_symb_dB = EbN0_dB(idx) + 10*log10(k);
        y = awgn(symbols .* H, snr_symb_dB, 'measured');

        % Ecualización con CSI perfecto
        y_eq = y ./ H;

        % Detección
        decoded_bits = zeros(num_bits, 1);
        for n = 1:num_symbols
            distances = abs(y_eq(n) - mapping).^2;
            [~, idx_min] = min(distances);
            decoded_bits(2*n-1:2*n) = bit_map(idx_min, :);
        end

        ber_total(idx) = ber_total(idx) + sum(decoded_bits' ~= bits);
    end
end

% BER promedio
ber_avg = ber_total / (num_runs * num_bits);

% BER teórica Rayleigh
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));

% Gráfica BER
figure;
semilogy(EbN0_dB, ber_rayleigh_theory, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_dB, ber_avg, 'r*-', 'LineWidth', 2);
grid on;
legend('Teórica Rayleigh', 'Sim. Multipath + AWGN (CSI perfecto)');
xlabel('E_b/N_0 [dB]'); ylabel('BER');
title('BER para QPSK en canal Multipath + AWGN');

% --- Visualización de constelaciones ---
% Nuevos bits
bits = randi([0 1], 1, num_bits);
bit_pairs = reshape(bits, 2, []).';
[~, indices] = ismember(bit_pairs, bit_map, 'rows');
symbols = mapping(indices).';

% Canal
t = linspace(0, 1, num_symbols);
an = ones(1,L)/sqrt(L);
thetan = 2*pi*rand(1,L);
fDn = fd_max * cos(2*pi*rand(1,L));
H = zeros(1,num_symbols);
for l = 1:L
    H = H + an(l)*exp(1j*(thetan(l) - 2*pi*fDn(l)*t));
end
H(abs(H) < 1e-3) = 1e-3;

% Ruido y señal
snr_example_dB = 30;
y_rx = awgn(symbols .* H, snr_example_dB + 10*log10(k), 'measured');
y_eq = y_rx ./ H;

% Gráficas de constelación
figure;
subplot(1,3,1); plot(real(symbols), imag(symbols), 'bo');
title('Constelación Transmitida'); axis equal; grid on;
xlabel('Re'); ylabel('Im');

subplot(1,3,2); plot(real(y_rx), imag(y_rx), 'rx');
title('Recibida sin ecualizar'); axis equal; grid on;
xlabel('Re'); ylabel('Im');

subplot(1,3,3); plot(real(y_eq), imag(y_eq), 'gx');
title('Recibida tras ecualización (CSI perfecto)');
axis equal; grid on; xlabel('Re'); ylabel('Im');