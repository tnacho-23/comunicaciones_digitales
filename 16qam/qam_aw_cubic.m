clc; clear all; close all;

% Parámetros fijos
num_bits = 1e5;
k = 4;  % 16QAM usa 4 bits por símbolo
EbN0_dB = -2:2:30;
num_runs = 21;

% Mapeo 16QAM Gray
re = [-3 -1 3 1];  % orden para código Gray
im = [-3 -1 3 1];
[re_grid, im_grid] = meshgrid(re, im);
mapping = (re_grid(:) + 1j*im_grid(:)) / sqrt(10);  % Normalización de energía

% Bit map Gray 16QAM (orden compatible con re/im)
bit_map = [0 0 0 0; 0 0 0 1; 0 0 1 1; 0 0 1 0;
           0 1 1 0; 0 1 1 1; 0 1 0 1; 0 1 0 0;
           1 1 0 0; 1 1 0 1; 1 1 1 1; 1 1 1 0;
           1 0 1 0; 1 0 1 1; 1 0 0 1; 1 0 0 0];

snr_plot_vals = [30];
L_vals = [10];
v_kmh_vals = [50];
fc_vals = [700e6, 3.5e9];
N_vals = [4, 9];
pilot_sym = 1 + 1j;
colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y', 'c'];
styles = {'-', '--', '-.', ':'};

for iN = 1:length(N_vals)
    N = N_vals(iN);
    saved_constellations = struct();
    figure;
    legend_entries = {};
    plot_idx = 1;

    for use_hamming = [false true]
        if use_hamming
            tag = '(Hamming)';
        else
            tag = '(sin codificación)';
        end

        for iL = 1:length(L_vals)
            for iv = 1:length(v_kmh_vals)
                for ifc = 1:length(fc_vals)
                    L = L_vals(iL);
                    v_kmh = v_kmh_vals(iv);
                    fc = fc_vals(ifc);
                    v = v_kmh / 3.6;
                    lambda = 3e8 / fc;
                    fd_max = v / lambda;

                    if use_hamming
                        rem = mod(num_bits, 4);
                        bits_padded = num_bits + (4 - rem)*(rem ~= 0);
                        num_bits_encoded = bits_padded * 7 / 4;
                    else
                        bits_padded = num_bits;
                        num_bits_encoded = num_bits;
                    end

                    num_symbols = num_bits_encoded / k;
                    ber_total = zeros(1, length(EbN0_dB));

                    for idx = 1:length(EbN0_dB)
                        EbN0_dB_val = EbN0_dB(idx);
                        EbN0 = 10^(EbN0_dB_val/10);
                        N0 = 1/(2*k*EbN0);

                        for run = 1:num_runs
                            bits_orig = randi([0 1], 1, num_bits);
                            if use_hamming
                                rem = mod(length(bits_orig), 4);
                                if rem ~= 0
                                    bits_orig = [bits_orig, zeros(1, 4 - rem)];
                                end
                                bits = hamming74_encode(bits_orig);
                            else
                                bits = bits_orig;
                            end

                            bit_quads = reshape(bits, k, []).';
                            indices = zeros(size(bit_quads, 1), 1);
                            for b = 1:size(bit_map, 1)
                                mask = all(bit_quads == bit_map(b,:), 2);
                                indices(mask) = b;
                            end
                            data_symbols = mapping(indices).';

                            tx_symbols = [];
                            pilot_pos = [];
                            data_pos = [];
                            i = 1;
                            while i <= length(data_symbols)
                                tx_symbols = [tx_symbols, pilot_sym];
                                pilot_pos = [pilot_pos, length(tx_symbols)];
                                if i+N-1 <= length(data_symbols)
                                    tx_symbols = [tx_symbols, data_symbols(i:i+N-1)];
                                    data_pos = [data_pos, length(tx_symbols)-N+1 : length(tx_symbols)];
                                    i = i + N;
                                else
                                    remaining = length(data_symbols) - i + 1;
                                    tx_symbols = [tx_symbols, data_symbols(i:end)];
                                    data_pos = [data_pos, length(tx_symbols)-remaining+1 : length(tx_symbols)];
                                    break;
                                end
                            end

                            len = length(tx_symbols);
                            t = linspace(0, 1, len);
                            an = ones(1,L)/sqrt(L);
                            thetan = 2*pi*rand(1,L);
                            fDn = fd_max * cos(2*pi*rand(1,L));
                            H_full = zeros(1,len);
                            for l = 1:L
                                H_full = H_full + an(l)*exp(1j*(thetan(l) - 2*pi*fDn(l)*t));
                            end
                            H_full(abs(H_full) < 1e-3) = 1e-3;
                            noise = sqrt(N0)*(randn(1,len) + 1j*randn(1,len));
                            y = tx_symbols .* H_full + noise;

                            rx_pilots = y(pilot_pos);
                            csi_est = rx_pilots / pilot_sym;
                            csi_interp_full = interp1(pilot_pos, csi_est, data_pos, 'spline');
                            rx_data = y(data_pos);
                            y_eq = rx_data ./ csi_interp_full;

                            decoded_bits = zeros(length(bits),1);
                            for n = 1:length(data_symbols)
                                distances = abs(y_eq(n) - mapping).^2;
                                [~, idx_min] = min(distances);
                                decoded_bits(k*n-k+1:k*n) = bit_map(idx_min,:).';
                            end

                            if use_hamming
                                rem = mod(length(decoded_bits), 7);
                                if rem ~= 0
                                    decoded_bits = [decoded_bits; zeros(7 - rem, 1)];
                                end
                                bits_final = hamming74_decode(decoded_bits.');
                                bits_final = bits_final(1:num_bits);
                                ber_total(idx) = ber_total(idx) + sum(bits_final ~= bits_orig);
                            else
                                ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);
                            end

                            if run == 1 && iL == 1 && iv == 1 && ifc == 1 && ismember(EbN0_dB_val, snr_plot_vals)
                                if EbN0_dB_val < 0
                                    s = sprintf('snr_m%d', abs(EbN0_dB_val));
                                else
                                    s = sprintf('snr_%d', EbN0_dB_val);
                                end
                                saved_constellations.(s).tx = data_symbols;
                                saved_constellations.(s).rx_awgn = data_symbols + noise(1:length(data_symbols));
                                saved_constellations.(s).rx_multipath = rx_data;
                                saved_constellations.(s).rx_eq = y_eq;
                            end
                        end
                    end

                    ber_avg = ber_total / (num_runs * num_bits);
                    style = styles{mod(plot_idx-1,length(styles))+1};
                    color = colors(mod(plot_idx-1,length(colors))+1);
                    semilogy(EbN0_dB, ber_avg, [color style], 'LineWidth', 1.8); hold on;
                    legend_entries{end+1} = sprintf('L=%d, v=%dkm/h, fc=%.1fGHz %s', L, v_kmh, fc/1e9, tag);
                    plot_idx = plot_idx + 1;
                end
            end
        end
    end

    grid on;
    legend(legend_entries, 'Location', 'southwest');
    xlabel('E_b/N_0 [dB]');
    ylabel('BER');
    title(sprintf('BER para 16QAM con N = %d símbolos entre pilotos (CSI por pilotos)', N+1));
end

% Funciones Hamming(7,4)
function encoded = hamming74_encode(bits)
    G = [1 0 0 0 1 1 0;
         0 1 0 0 1 0 1;
         0 0 1 0 1 0 0;
         0 0 0 1 0 1 1];
    k = 4;
    num_blocks = length(bits)/k;
    bits_reshaped = reshape(bits, k, num_blocks).';
    encoded = mod(bits_reshaped * G, 2);
    encoded = reshape(encoded.', 1, []);
end

function decoded = hamming74_decode(bits)
    H = [1 1 1 0 1 0 0;
         1 0 0 1 0 1 0;
         0 1 0 1 0 0 1];
    n = 7;
    num_blocks = length(bits)/n;
    bits_reshaped = reshape(bits, n, num_blocks).';
    decoded = zeros(num_blocks, 4);
    for i = 1:num_blocks
        r = bits_reshaped(i,:);
        syndrome = mod(H * r.', 2);
        syndrome_dec = bi2de(syndrome.', 'left-msb');
        if syndrome_dec ~= 0 && syndrome_dec <= n
            r(syndrome_dec) = mod(r(syndrome_dec) + 1, 2);
        end
        decoded(i,:) = r(1:4);
    end
    decoded = reshape(decoded.', 1, []);
end
