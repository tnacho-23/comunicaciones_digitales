clc; clear all; close all;

% Parámetros
num_bits = 1e5;
k = 2; % bits por símbolo para QPSK
num_data_symbols = num_bits / k;
pilot_spacing = 4;

% Estimar número total de símbolos con pilotos
num_symbols = ceil(num_data_symbols / (1 - 1/pilot_spacing));

% Recalcular número de pilotos y datos
pilot_positions = 1:pilot_spacing:num_symbols;
num_pilots = length(pilot_positions);
data_positions = setdiff(1:num_symbols, pilot_positions);
num_data_symbols = length(data_positions); % ajustar datos reales

% Ajustar número de bits en base a los símbolos de datos reales
num_bits = num_data_symbols * k;

EbN0_dB = -2:1:30;
L = 50; % número de trayectorias
v_kmh = 120; fc = 700e6;
v = v_kmh / 3.6;
lambda = 3e8 / fc;
fd_max = v / lambda;
num_runs = 1;

% Mapeo QPSK Gray
bit_map = [0 0; 0 1; 1 1; 1 0]; % Codificación Gray
mapping = [1+1j; -1+1j; -1-1j; 1-1j] / sqrt(2); % Símbolos QPSK
pilot_symbol = 1 + 1j; % símbolo piloto conocido

ber_total = zeros(1, length(EbN0_dB));

for idx = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(idx)/10);
    for run = 1:num_runs
        % Generar bits aleatorios
        bits = randi([0 1], 1, num_bits);
        bit_pairs = reshape(bits, 2, []).';
        [~, indices] = ismember(bit_pairs, bit_map, 'rows');
        data_symbols = mapping(indices).';

        % Insertar pilotos
        symbols = zeros(1, num_symbols);
        symbols(pilot_positions) = pilot_symbol;
        symbols(data_positions) = data_symbols;

        % Canal Rayleigh con L trayectorias
        t = linspace(0, 1, num_symbols);
        an = ones(1,L)/sqrt(L);
        thetan = 2*pi*rand(1,L);
        fDn = fd_max * cos(2*pi*rand(1,L));
        H = zeros(1,num_symbols);
        for l = 1:L
            H = H + an(l)*exp(1j*(thetan(l) - 2*pi*fDn(l)*t));
        end
        H(abs(H) < 1e-3) = 1e-3;

        % Transmisión con ruido
        snr_symb_dB = EbN0_dB(idx) + 10*log10(k);
        y = awgn(symbols .* H, snr_symb_dB, 'measured');

        % Estimación del canal en pilotos
        H_est = zeros(1,num_symbols);
        H_est(pilot_positions) = y(pilot_positions) ./ pilot_symbol;

        % Interpolación cúbica para el canal en símbolos de datos
        H_est(data_positions) = interp1(pilot_positions, H_est(pilot_positions), data_positions, 'pchip', 'extrap');

        % Ecualización
        y_eq = y ./ H_est;

        % Demodulación de los símbolos de datos
        decoded_bits = zeros(num_bits, 1);
        for n = 1:length(data_positions)
            symb = y_eq(data_positions(n));
            distances = abs(symb - mapping).^2;
            [~, idx_min] = min(distances);
            decoded_bits(2*n-1:2*n) = bit_map(idx_min, :);
        end

        % Cálculo de BER
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
legend('Teórica Rayleigh', 'Sim. Estimación cúbica con pilotos');
xlabel('E_b/N_0 [dB]'); ylabel('BER');
title('BER para QPSK en canal Multipath + AWGN con interpolación cúbica');

%% --- Visualización de constelaciones ---
% Nuevos bits
bits = randi([0 1], 1, num_bits);
bit_pairs = reshape(bits, 2, []).';
[~, indices] = ismember(bit_pairs, bit_map, 'rows');
symbols_data = mapping(indices).';

% Insertar pilotos
symbols = zeros(1, num_symbols);
symbols(pilot_positions) = pilot_symbol;
symbols(data_positions) = symbols_data;

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

% Gráficas de constelación (solo símbolos de datos)
figure;
subplot(1,3,1); plot(real(symbols(data_positions)), imag(symbols(data_positions)), 'bo');
title('Constelación Transmitida'); axis equal; grid on;
xlabel('Re'); ylabel('Im');

subplot(1,3,2); plot(real(y_rx(data_positions)), imag(y_rx(data_positions)), 'rx');
title('Recibida sin ecualizar'); axis equal; grid on;
xlabel('Re'); ylabel('Im');

subplot(1,3,3); plot(real(y_eq(data_positions)), imag(y_eq(data_positions)), 'gx');
title('Recibida tras ecualización (CSI perfecto)');
axis equal; grid on; xlabel('Re'); ylabel('Im');
