clc; clear; close all;

% Parámetros fijos
num_bits = 1e6;
k = 2;
num_symbols = num_bits / k;
EbN0_dB = -2:2:30;
EbN0_lin = 10.^(EbN0_dB / 10);
num_runs = 10;
pilot_interval = 5;
mapping = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 0; 1 1];
pilot_symbol = 1 + 1j;

% Variables a recorrer
L_vals = [5];
v_kmh_vals = [30];
fc_vals = [700e6];

% Gráfica
colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y', 'c'];
styles = {'-', '--', '-.', ':'};
figure;
legend_entries = {};
plot_idx = 1;

for iL = 1:length(L_vals)
    for iv = 1:length(v_kmh_vals)
        for ifc = 1:length(fc_vals)
            L = L_vals(iL);
            v_kmh = v_kmh_vals(iv);
            fc = fc_vals(ifc);

            v = v_kmh / 3.6;
            lambda = 3e8 / fc;
            fd_max = v / lambda;
            an = ones(1, L) / sqrt(L);

            num_pilots = ceil(num_symbols / (pilot_interval - 1));
            total_symbols = num_symbols + num_pilots;
            pilot_indices = pilot_interval:pilot_interval:total_symbols;
            data_indices = setdiff(1:total_symbols, pilot_indices);
            data_indices = data_indices(1:num_symbols); % asegurar tamaño correcto

            t = linspace(0, 1, total_symbols);
            ber_total = zeros(1, length(EbN0_dB));

            for run = 1:num_runs
                bits = randi([0 1], num_bits, 1);
                bit_pairs = reshape(bits, 2, []).';
                indices = bit_pairs(:, 1) * 2 + bit_pairs(:, 2) + 1;
                symbols = mapping(indices).';

                symbols_tx = zeros(1, total_symbols);
                symbols_tx(pilot_indices) = pilot_symbol;
                symbols_tx(data_indices) = symbols;

                thetan = 2 * pi * rand(1, L);
                fDn = fd_max * cos(2 * pi * rand(1, L));
                H = sum(an.' .* exp(1j * (thetan.' - 2 * pi * fDn.' * t)), 1);
                H(abs(H) < 1e-3) = 1e-3;

                for idx_ebn0 = 1:length(EbN0_dB)
                    N0 = 1 / (2 * k * EbN0_lin(idx_ebn0));
                    noise = sqrt(N0) * (randn(1, total_symbols) + 1j * randn(1, total_symbols));
                    y = symbols_tx .* H + noise;

                    % Estimación de canal con interpolación FFT
                    H_est = zeros(1, total_symbols);
                    H_est(pilot_indices) = y(pilot_indices) / pilot_symbol;

                    % Suavizar la transición (relleno suave)
                    if all(H_est == 0)
                        warning('Estimación H_est vacía. Problema en pilotos.');
                        continue;
                    end

                    H_est = ifft(fft(H_est));
                    H_est(abs(H_est) < 1e-3) = 1e-3;

                    % Ecualización ZF
                    y_eq = y ./ H_est;

                    % Decodificación
                    decoded_bits = zeros(num_bits, 1);
                    for n = 1:num_symbols
                        idx = data_indices(n);
                        [~, idx_min] = min(abs(y_eq(idx) - mapping).^2);
                        decoded_bits(2 * n - 1 : 2 * n) = bit_map(idx_min, :).';
                    end

                    if length(decoded_bits) ~= length(bits)
                        error('Desalineación entre decoded_bits y bits');
                    end

                    ber_total(idx_ebn0) = ber_total(idx_ebn0) + sum(decoded_bits ~= bits);
                end
            end

            ber_avg = ber_total / (num_runs * num_bits);

            % Gráfico
            style = styles{mod(plot_idx - 1, length(styles)) + 1};
            color = colors(mod(plot_idx - 1, length(colors)) + 1);
            semilogy(EbN0_dB, ber_avg, [color style], 'LineWidth', 1.8); hold on;
            legend_entries{end + 1} = sprintf('L=%d, v=%dkm/h, fc=%.1fGHz', L, v_kmh, fc / 1e9);
            plot_idx = plot_idx + 1;
        end
    end
end

% Curva teórica Rayleigh
ber_rayleigh_theory = 0.5 * (1 - sqrt(EbN0_lin ./ (EbN0_lin + 1)));
semilogy(EbN0_dB, ber_rayleigh_theory, 'k--', 'LineWidth', 2);
legend_entries{end + 1} = 'Teórica Rayleigh';

grid on;
legend(legend_entries, 'Location', 'southwest');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
title('BER QPSK con estimación FFT y ZF');
