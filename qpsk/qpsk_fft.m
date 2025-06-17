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

% Pilotos
pilot_symbol = 1 + 1j;
pilot_interval = 5; 

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

        % Insertar pilotos cada 'pilot_interval' símbolos
        total_len = num_symbols + ceil(num_symbols/pilot_interval);
        tx_frame = zeros(1, total_len);
        is_pilot = false(1, total_len);
        data_ptr = 1;

        for i = 1:total_len
            if mod(i-1, pilot_interval + 1) == 0
                tx_frame(i) = pilot_symbol;
                is_pilot(i) = true;
            else
                tx_frame(i) = symbols(data_ptr);
                data_ptr = data_ptr + 1;
            end
        end

        % Canal
        t = linspace(0, 1, length(tx_frame));
        an = ones(1,L)/sqrt(L);
        thetan = 2*pi*rand(1,L);
        fDn = fd_max * cos(2*pi*rand(1,L));
        H = zeros(1,length(tx_frame));
        for l = 1:L
            H = H + an(l)*exp(1j*(thetan(l) - 2*pi*fDn(l)*t));
        end

        noise = sqrt(N0)*(randn(1,length(tx_frame)) + 1j*randn(1,length(tx_frame)));
        y = tx_frame .* H + noise;

        % Estimación del canal en pilotos
        H_est = zeros(1, length(tx_frame));
        pilot_pos = find(is_pilot);
        H_est(pilot_pos) = y(pilot_pos) / pilot_symbol;

        % Interpolación
        data_pos = find(~is_pilot);
        H_est(data_pos) = interp1(pilot_pos, H_est(pilot_pos), data_pos, 'linear', 'extrap');

        % Ecualización
        y_eq_all = y ./ H_est;

        % Extraer datos
        y_data = y_eq_all(~is_pilot);

        % Detección
        decoded_bits = zeros(num_bits, 1);
        for n = 1:num_symbols
            distances = abs(y_data(n) - mapping).^2;
            [~, idx_min] = min(distances);
            decoded_bits(2*n-1:2*n) = bit_map(idx_min, :);
        end

        ber_total(idx) = ber_total(idx) + sum(decoded_bits' ~= bits);

        % ===== Mostrar constelaciones si Eb/N0 = 30 dB (último idx) =====
        if EbN0_dB(idx) == 30 && run == 1
            figure;
            subplot(1,3,1); plot(real(symbols), imag(symbols), 'bo');
            title('Constelación Transmitida');
            xlabel('Re'); ylabel('Im'); axis equal; grid on;

            y_rx_data = y(~is_pilot);
            subplot(1,3,2); plot(real(y_rx_data), imag(y_rx_data), 'rx');
            title('Recibida sin ecualizar');
            xlabel('Re'); ylabel('Im'); axis equal; grid on;

            subplot(1,3,3); plot(real(y_data), imag(y_data), 'gx');
            title('Ecualizada (est. canal)');
            xlabel('Re'); ylabel('Im'); axis equal; grid on;
        end
    end
end

% Resultados
ber_avg = ber_total / (num_runs * num_bits);
% disp(ber_avg)

% Teórica Rayleigh
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));

% Gráfica de BER
figure;
semilogy(EbN0_dB, ber_rayleigh_theory, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_dB, ber_avg, 'ro-', 'LineWidth', 2);
grid on;
legend('Teórica Rayleigh', ['Sim. con piloto cada ', num2str(pilot_interval+1), ' símbolos']);
xlabel('E_b/N_0 [dB]'); ylabel('BER');
title('BER para QPSK con estimación de canal por pilotos');
