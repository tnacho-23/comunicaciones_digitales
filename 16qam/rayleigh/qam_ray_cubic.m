clc; clear all; close all;

% Distancias entre pilotos a evaluar
N_vals = [1, 4, 9, 19];

for iN = 1:length(N_vals)
    N = N_vals(iN);
    fprintf('\n--- Simulación con piloto cada %d símbolos (N = %d) ---\n', N, N);

    % Parámetros fijos
    num_bits = 1e3;
    k = 4;  % 16-QAM usa 4 bits por símbolo
    num_symbols = num_bits / k;
    EbN0_dB = -2:2:30;
    num_runs = 10;

    % Constelación 16-QAM normalizada
    re_vals = [-3, -1, 1, 3];
    im_vals = [-3, -1, 1, 3];
    [re_grid, im_grid] = meshgrid(re_vals, im_vals);
    mapping = (re_grid(:) + 1j*im_grid(:)) / sqrt(10);
    bit_map = de2bi(0:15, 4, 'left-msb');

    snr_plot_vals = [-2, 0, 10, 30];
    saved_constellations = struct();

    % Canal plano: valores arbitrarios para L, v, fc (no se usarán)
    L_vals = [5];
    v_kmh_vals = [0];
    fc_vals = [2.4e9];

    pilot_sym = mapping(1);

    colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y', 'c'];
    styles = {'-', '--', '-.', ':'};

    figure;
    legend_entries = {};
    plot_idx = 1;

    for iL = 1:length(L_vals)
        for iv = 1:length(v_kmh_vals)
            for ifc = 1:length(fc_vals)
                L = L_vals(iL);
                ber_total = zeros(1, length(EbN0_dB));

                for idx = 1:length(EbN0_dB)
                    EbN0_dB_val = EbN0_dB(idx);
                    EbN0 = 10^(EbN0_dB_val/10);
                    N0 = 1/(2*k*EbN0);

                    for run = 1:num_runs
                        bits = randi([0 1], 1, num_bits);
                        bit_quads = reshape(bits, 4, []).';
                        indices = bi2de(bit_quads, 'left-msb') + 1;
                        data_symbols = mapping(indices).';

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

                        % Canal Rayleigh plano (independiente por símbolo)
                        H_full = (randn(1,len) + 1j*randn(1,len)) / sqrt(2);
                        H_full(abs(H_full) < 1e-3) = 1e-3;

                        noise = sqrt(N0)*(randn(1,len) + 1j*randn(1,len));
                        y = tx_symbols .* H_full + noise;

                        % Estimación CSI (spline)
                        rx_pilots = y(pilot_pos);
                        csi_est = rx_pilots / pilot_sym;
                        csi_interp_full = interp1(pilot_pos, csi_est, data_pos, 'spline');

                        rx_data = y(data_pos);
                        y_eq = rx_data ./ csi_interp_full;

                        decoded_bits = zeros(num_bits,1);
                        for n = 1:num_symbols
                            distances = abs(y_eq(n) - mapping).^2;
                            [~, idx_min] = min(distances);
                            decoded_bits(4*n-3:4*n) = bit_map(idx_min,:).';
                        end
                        ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);

                        should_save = (run == 1) && (iL == 1) && (iv == 1) && (ifc == 1);
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
                            fprintf('Guardando constelaciones para %s (N = %d)\n', s, N);
                        end
                    end
                end

                ber_avg = ber_total / (num_runs * num_bits);

                style = styles{mod(plot_idx-1,length(styles))+1};
                color = colors(mod(plot_idx-1,length(colors))+1);
                semilogy(EbN0_dB, ber_avg, [color style], 'LineWidth', 1.8); hold on;

                legend_entries{end+1} = sprintf('Canal plano (L=%d)', L);
                plot_idx = plot_idx + 1;
            end
        end
    end

    % Curva teórica para 16-QAM en Rayleigh
    EbN0_lin = 10.^(EbN0_dB/10);
    M = 16;
    k = log2(M);
    gamma = (3/(M-1)) * EbN0_lin;
    Pb_16qam_rayleigh = (4/k) * (1 - 1/sqrt(M)) * 0.5 .* (1 - sqrt(gamma ./ (gamma + 1)));
    semilogy(EbN0_dB, Pb_16qam_rayleigh, 'k--', 'LineWidth', 2);
    legend_entries{end+1} = 'Teórica Rayleigh 16QAM';

    grid on;
    legend(legend_entries, 'Location', 'southwest');
    xlabel('E_b/N_0 [dB]');
    ylabel('BER');
    title(sprintf('BER para 16QAM (spline, N = %d, canal plano)', N));

    %% Constelaciones
    figure;
    snrs = snr_plot_vals;
    titles = {'1. Original', '2. Solo AWGN', '3. Rayleigh + AWGN', '4. Ecualizado'};

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
            fprintf('Advertencia: no se guardaron constelaciones para SNR = %d dB (N = %d)\n', snrs(i), N);
        end
    end

    sgtitle(sprintf('Constelaciones spline para SNR = -2, 0, 10, 30 dB (N = %d)', N));
end
