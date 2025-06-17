clc; clear all; close all;

% Parámetros
num_bits = 1000;
k = 2; % bits por símbolo para QPSK
num_symbols = num_bits / k;
EbN0_dB = -2:1:30;
L = 1; % número de reflexiones
v_kmh = 0; fc = 700e6;
v = v_kmh / 3.6;
lambda = 3e8 / fc;
fd_max = v / lambda;
num_runs = 1;

% Mapeo QPSK (Gray)
mapping = [1+1j; -1+1j; -1-1j; 1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 1; 1 0];

ber_total = zeros(1, length(EbN0_dB));

for idx = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(idx)/10);
    N0 = 1/(2*k*EbN0);

    for run = 1:num_runs
        % Transmisión
        bits = randi([0 1], 1, num_bits);
        bit_pairs = reshape(bits, 2, []).';
        indices = bit_pairs(:,1)*2 + bit_pairs(:,2) + 1;
        symbols = mapping(indices).';

        t = linspace(0, 1, num_symbols);
        an = ones(1,L)/sqrt(L);
        thetan = 2*pi*rand(1,L);
        fDn = fd_max * cos(2*pi*rand(1,L));
        H = zeros(1,num_symbols);
        for l = 1:L
            H = H + an(l)*exp(1j*(thetan(l) - 2*pi*fDn(l)*t));
        end

        % Evitar divisiones por valores cercanos a cero
        H(abs(H) < 1e-3) = 1e-3;

        noise = sqrt(N0)*(randn(1,num_symbols) + 1j*randn(1,num_symbols));
        y = symbols .* H + noise;

        % Ecualización (CSI perfecto)
        y_eq = y ./ H;

        %figure;
        %plot(real(y_eq), imag(y_eq), 'gx');
        %title(['Constelación ecualizada - E_b/N_0 = ', num2str(EbN0_dB(idx)), ' dB']);
        %xlabel('Re'); ylabel('Im');
        %axis equal; grid on;


        % Detección
        decoded_bits = zeros(num_bits, 1);
        for n = 1:num_symbols
            distances = abs(y_eq(n) - mapping).^2; % Distancias al cuadrado
            [~, idx_min] = min(distances);         % Índice del símbolo más cercano
            decoded_bits(2*n-1:2*n) = bit_map(idx_min, :); % Recuperar los bits originales
        end


        ber_total(idx) = ber_total(idx) + sum(decoded_bits' ~= bits);
    end
end

%disp(ber_total)
% Normalización
ber_avg = ber_total / (num_runs * num_bits);
%disp(ber_avg)

% Teórica Rayleigh
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));

% Gráfica de BER
figure;
semilogy(EbN0_dB, ber_rayleigh_theory, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_dB, ber_avg, 'ro-', 'LineWidth', 2);
grid on;
legend('Teórica Rayleigh', 'Sim. Multipath + AWGN (CSI perfecto)');
xlabel('E_b/N_0 [dB]'); ylabel('BER');
title('BER para QPSK en canal Multipath + AWGN');

% Constelaciones
bits = randi([0 1], 1, num_bits);
bit_pairs = reshape(bits, 2, []).';
indices = bit_pairs(:,1)*2 + bit_pairs(:,2) + 1;
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

% Evitar divisiones por cero en H
H(abs(H) < 1e-3) = 1e-3;

% Ruido y señal
snr_example_dB = 30;
N0 = 1/(2*k*10^(snr_example_dB/10));
noise = sqrt(N0)*(randn(1,num_symbols) + 1j*randn(1,num_symbols));
y_rx = symbols .* H + noise;
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
