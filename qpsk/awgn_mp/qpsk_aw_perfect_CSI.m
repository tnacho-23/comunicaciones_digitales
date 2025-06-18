clc; clear all; close all;

% Parámetros
num_bits = 1e5;     % Más bits para buena estadística
k = 2;              % Bits por símbolo (QPSK)
num_symbols = num_bits / k;
EbN0_dB = -2:1:30;
L = 50;              % Número de reflexiones
v_kmh = 120;         % Velocidad móvil
fc = 700e6;
v = v_kmh / 3.6;
lambda = 3e8 / fc;
fd_max = v / lambda;
num_runs = 10;

% Mapeo QPSK (Gray correcto)
mapping = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 0; 1 1];  % Bits que corresponden a cada símbolo

ber_total = zeros(1, length(EbN0_dB));

for idx = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(idx)/10);
    N0 = 1/(2*k*EbN0);

    for run = 1:num_runs
        bits = randi([0 1], 1, num_bits);
        bit_pairs = reshape(bits, 2, []).';

        indices = bit_pairs(:,1)*2 + bit_pairs(:,2) + 1; % b1 MSB, b2 LSB
        symbols = mapping(indices).';

        % Canal Rayleigh multipath
        t = linspace(0, 1, num_symbols);
        an = ones(1,L)/sqrt(L);
        thetan = 2*pi*rand(1,L);
        fDn = fd_max * cos(2*pi*rand(1,L));
        H = zeros(1,num_symbols);
        for l = 1:L
            H = H + an(l)*exp(1j*(thetan(l) - 2*pi*fDn(l)*t));
        end
        H(abs(H) < 1e-3) = 1e-3;

        noise = sqrt(N0)*(randn(1,num_symbols) + 1j*randn(1,num_symbols));
        y = symbols .* H + noise;
        y_eq = y ./ H;

        decoded_bits = zeros(num_bits,1);
        for n = 1:num_symbols
            distances = abs(y_eq(n) - mapping).^2;
            [~, idx_min] = min(distances);
            decoded_bits(2*n-1:2*n) = bit_map(idx_min,:).';
        end

        ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);
    end
end

ber_avg = ber_total / (num_runs * num_bits);

% BER teórica
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));

% Gráfica BER
figure;
semilogy(EbN0_dB, ber_rayleigh_theory, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_dB, ber_avg, 'ro-', 'LineWidth', 2);
grid on;
legend('Teórica Rayleigh', 'Simulación Multipath + AWGN (CSI perfecto)');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
title('BER para QPSK en canal Multipath + AWGN');

% Constelaciones ejemplo a 30 dB
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
H(abs(H) < 1e-3) = 1e-3;

snr_example_dB = 30;
N0 = 1/(2*k*10^(snr_example_dB/10));
noise = sqrt(N0)*(randn(1,num_symbols) + 1j*randn(1,num_symbols));
y_rx = symbols .* H + noise;
y_eq = y_rx ./ H;

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
