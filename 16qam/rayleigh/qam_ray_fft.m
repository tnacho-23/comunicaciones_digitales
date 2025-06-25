clc; clear all; close all;

% Parámetros fijos
num_bits = 1e5; % Más bits para buena estadística
k = 4; % 16-QAM -> 4 bits/símbolo
num_symbols = num_bits / k;
EbN0_dB = -2:2:30;
num_runs = 21;
snr_plot_vals = [-2, 0, 10, 30];
fc_vals = [700e6, 3.5e9];  % Se mantiene

% Pilotos a probar
N_vals = [1, 4, 9, 19];
pilot_sym = 1 + 1j;

% Colores y estilos
colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y', 'c'];
styles = {'-', '--', '-.', ':'};

% Mapeo 16-QAM con Gray coding (normalizado)
re = [-3 -1 +3 +1];
im = [-3 -1 +3 +1];
[X, Y] = meshgrid(re, im);
constellation = X(:) + 1j * Y(:);
constellation = constellation / sqrt(mean(abs(constellation).^2)); % Normalizar potencia

% Mapa de bits correspondiente (Gray code para 16-QAM)
bit_map = de2bi([0:15], 4, 'left-msb'); % 16-QAM: 4 bits por símbolo

for iN = 1:length(N_vals)
    N = N_vals(iN);
    saved_constellations = struct();
    figure;
    legend_entries = {};
    plot_idx = 1;

    for ifc = 1:length(fc_vals)
        fc = fc_vals(ifc);

        ber_total = zeros(1, length(EbN0_dB));

        for idx = 1:length(EbN0_dB)
            EbN0_dB_val = EbN0_dB(idx);
            EbN0 = 10^(EbN0_dB_val/10);
            N0 = 1/(2*k*EbN0);

            for run = 1:num_runs
                % Bits aleatorios y mapeo 16-QAM
                bits = randi([0 1], 1, num_bits);
                bit_groups = reshape(bits, k, []).';
                dec_vals = bi2de(bit_groups, 'left-msb') + 1;
                data_symbols = constellation(dec_vals).';

                % Insertar pilotos
                tx_symbols = [];
                pilot_pos = [];
                data_pos = [];
                i = 1;
                while i <= num_symbols
                    tx_symbols = [tx_symbols, pilot_sym];
                    pilot_pos = [pilot_pos, length(tx_symbols)];

                    if i+N-1 <= num_symbols
                        tx_symbols = [tx_symbols, data_symbols(i:i+N-1)];
                        data_pos = [data_pos, length(tx_symbols)-N+1 : length(tx_symbols)];
                        i = i + N;
                    else
                        remaining = num_symbols - i + 1;
                        tx_symbols = [tx_symbols, data_symbols(i:end)];
                        data_pos = [data_pos, length(tx_symbols)-remaining+1 : length(tx_symbols)];
                        break;
                    end
                end

                len = length(tx_symbols);

                % Canal Rayleigh plano
                H_full = (randn(1,len) + 1j*randn(1,len)) / sqrt(2);
                H_full(abs(H_full) < 1e-3) = 1e-3;

                noise = sqrt(N0)*(randn(1,len) + 1j*randn(1,len));
                y = tx_symbols .* H_full + noise;

                % Estimación CSI por pilotos
                rx_pilots = y(pilot_pos);
                csi_est = rx_pilots / pilot_sym;
                csi_interp_full = interpft(csi_est, length(tx_symbols));

                % Ecualización
                rx_data = y(data_pos);
                csi_used = csi_interp_full(data_pos);
                y_eq = rx_data ./ csi_used;

                % Decodificación por distancia mínima
                decoded_bits = zeros(num_bits,1);
                for n = 1:num_symbols
                    distances = abs(y_eq(n) - constellation).^2;
                    [~, idx_min] = min(distances);
                    decoded_bits((n-1)*k+1 : n*k) = bit_map(idx_min, :).';
                end
                ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);

                % Guardar constelaciones
                should_save = (run == 1) && (ifc == 1);
                if should_save && ismember(EbN0_dB_val, snr_plot_vals)
                    if EbN0_dB_val < 0
                        s = sprintf('snr_m%d', abs(EbN0_dB_val));
                    else
                        s = sprintf('snr_%d', EbN0_dB_val);
                    end
                    saved_constellations.(s).tx = data_symbols;
                    saved_constellations.(s).rx_awgn = data_symbols + noise(1:num_symbols);
                    saved_constellations.(s).rx_multipath = rx_data;
                    saved_constellations.(s).rx_eq = y_eq;
                end
            end
        end

        ber_avg = ber_total / (num_runs * num_bits);
        style = styles{mod(plot_idx-1,length(styles))+1};
        color = colors(mod(plot_idx-1,length(colors))+1);
        semilogy(EbN0_dB, ber_avg, [color style], 'LineWidth', 1.8); hold on;

        legend_entries{end+1} = sprintf('fc = %.1f GHz', fc/1e9);
        plot_idx = plot_idx + 1;
    end

    % No hay curva teórica exacta simple para 16-QAM en Rayleigh
    EbN0_lin = 10.^(EbN0_dB/10);
    ber_rayleigh_theory = 3/8 * (1 - sqrt(EbN0_lin./(EbN0_lin+1))); % Aprox
    semilogy(EbN0_dB, ber_rayleigh_theory, 'k--', 'LineWidth', 2);
    legend_entries{end+1} = 'Aprox. Rayleigh';

    grid on;
    legend(legend_entries, 'Location', 'southwest');
    xlabel('E_b/N_0 [dB]');
    ylabel('BER');
    title(sprintf('BER para 16-QAM, pilotos cada %d datos (Rayleigh plano + AWGN)', N));

    %% Figura de constelaciones
    figure;
    snrs = snr_plot_vals;
    titles = {'1. Original', '2. Solo AWGN', '3. Canal + AWGN', '4. Ecualizado'};

    for i = 1:length(snrs)
        if snrs(i) < 0
            s = sprintf('snr_m%d', abs(snrs(i)));
        else
            s = sprintf('snr_%d', snrs(i));
        end

        if isfield(saved_constellations, s)
            data = saved_constellations.(s);

            subplot(4,4,(i-1)*4+1);
            plot(real(data.tx), imag(data.tx), 'bo'); axis equal; grid on;
            title(titles{1});
            ylabel(sprintf('SNR = %d dB', snrs(i)));

            subplot(4,4,(i-1)*4+2);
            plot(real(data.rx_awgn), imag(data.rx_awgn), 'ro'); axis equal; grid on;
            title(titles{2});

            subplot(4,4,(i-1)*4+3);
            plot(real(data.rx_multipath), imag(data.rx_multipath), 'mo'); axis equal; grid on;
            title(titles{3});

            subplot(4,4,(i-1)*4+4);
            plot(real(data.rx_eq), imag(data.rx_eq), 'go'); axis equal; grid on;
            title(titles{4});
        else
            fprintf('Advertencia: no se guardaron constelaciones para SNR = %d dB\n', snrs(i));
        end
    end

    sgtitle(sprintf('Constelaciones FFT (N = %d símbolos entre pilotos)', N));
end
