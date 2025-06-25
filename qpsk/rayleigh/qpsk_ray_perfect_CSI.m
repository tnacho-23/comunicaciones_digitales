clc; clear all; close all;

% Parámetros fijos
num_bits = 1e4;
k = 2;
num_symbols = num_bits / k;
EbN0_dB = -2:2:30;
num_runs = 21;
mapping = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 0; 1 1];

% Parámetros variables
fc_vals = [700e6, 3.5e9];

% Para graficar
colors = lines(length(fc_vals));
plot_idx = 1;
legend_entries = {};
constellation_plotted = false;

figure(1); % Figura para BER

for ifc = 1:length(fc_vals)

    fc = fc_vals(ifc);

    ber_total = zeros(1, length(EbN0_dB));

    for idx = 1:length(EbN0_dB)
        EbN0 = 10^(EbN0_dB(idx)/10);
        N0 = 1/(2*k*EbN0);

        for run = 1:num_runs
            bits = randi([0 1], 1, num_bits);
            bit_pairs = reshape(bits, 2, []).';
            indices = bit_pairs(:,1)*2 + bit_pairs(:,2) + 1;
            symbols = mapping(indices).';

            % Canal plano Rayleigh
            H = (randn(1,num_symbols) + 1j*randn(1,num_symbols))/sqrt(2);

            noise = sqrt(N0)*(randn(1,num_symbols) + 1j*randn(1,num_symbols));
            y = symbols .* H + noise;
            y_eq = y ./ H;

            % Graficar constelaciones (una sola vez)
            if ~constellation_plotted && ifc == 1 && run == 1 && idx == length(EbN0_dB)
                figure(2); % Figura para constelaciones
                sgtitle(sprintf('Constelaciones'));

                subplot(1,3,1);
                plot(real(symbols), imag(symbols), 'o');
                axis equal; grid on;
                title('Transmitido (ideal)');
                xlabel('Re'); ylabel('Im');

                subplot(1,3,2);
                plot(real(y), imag(y), 'x');
                axis equal; grid on;
                title('Recibido sin ecualizar');
                xlabel('Re'); ylabel('Im');

                subplot(1,3,3);
                plot(real(y_eq), imag(y_eq), '+');
                axis equal; grid on;
                title('Recibido tras ecualizar');
                xlabel('Re'); ylabel('Im');

                constellation_plotted = true;
            end

            % Decodificación
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

    % Graficar BER
    figure(1);
    semilogy(EbN0_dB, ber_avg, '-', 'LineWidth', 1.8, ...
             'Color', colors(plot_idx,:)); hold on;

    legend_entries{end+1} = sprintf('fc=%.1fGHz', fc/1e9);
    plot_idx = plot_idx + 1;
end

% BER teórica para canal Rayleigh plano
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));
figure(1);
semilogy(EbN0_dB, ber_rayleigh_theory, 'k--', 'LineWidth', 2);
legend_entries{end+1} = 'Teórica Rayleigh';

grid on;
legend(legend_entries, 'Location', 'southwest');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
title('BER para QPSK en canal plano Rayleigh + AWGN (CSI perfecto)');
