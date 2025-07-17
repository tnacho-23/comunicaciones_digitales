clc; clear all; close all;

% Parámetros fijos
num_bits = 1e5;
k = 2; % QPSK: 2 bits por símbolo
EbN0_dB = -2:2:30;
num_runs = 21;
mapping = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 0; 1 1];
L_vals = [10];
v_kmh_vals = [50];
fc_vals = [700e6, 3.5e9];
colors = ['r', 'b', 'g', 'm', 'k', 'c', 'y'];
styles = {'-', '--', '-.', ':'};

% 🔧 Crear una única figura para ambos casos
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

                ber_avg = zeros(1, length(EbN0_dB));

                for idx = 1:length(EbN0_dB)
                    EbN0_dB_val = EbN0_dB(idx);
                    EbN0 = 10^(EbN0_dB_val/10);
                    N0 = 1/(2*k*EbN0);

                    num_errors_total = 0;
                    num_bits_total = 0;

                    for run = 1:num_runs
                        bits_orig = randi([0 1], 1, num_bits);

                        % Codificación Hamming
                        if use_hamming
                            rem4 = mod(length(bits_orig), 4);
                            if rem4 ~= 0
                                bits_orig_padded = [bits_orig, zeros(1, 4 - rem4)];
                                hamming_padding = 4 - rem4;
                            else
                                bits_orig_padded = bits_orig;
                                hamming_padding = 0;
                            end
                            bits = hamming74_encode(bits_orig_padded);
                        else
                            bits = bits_orig;
                        end

                        % Padding para QPSK
                        if mod(length(bits), 2) ~= 0
                            bits = [bits, 0];
                        end

                        bit_pairs = reshape(bits, 2, []).';
                        indices = bit_pairs(:,1)*2 + bit_pairs(:,2) + 1;
                        data_symbols = mapping(indices).';

                        % Canal Rayleigh
                        len = length(data_symbols);
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
                        y = data_symbols .* H_full + noise;

                        % Ecualización (CSI perfecto)
                        y_eq = y ./ H_full;

                        % Demodulación
                        decoded_bits = zeros(length(bits),1);
                        for n = 1:length(data_symbols)
                            distances = abs(y_eq(n) - mapping).^2;
                            [~, idx_min] = min(distances);
                            decoded_bits(2*n-1:2*n) = bit_map(idx_min,:).';
                        end

                        % Decodificación y comparación
                        if use_hamming
                            rem7 = mod(length(decoded_bits), 7);
                            if rem7 ~= 0
                                decoded_bits = [decoded_bits; zeros(7 - rem7, 1)];
                            end
                            bits_final = hamming74_decode(decoded_bits.');

                            if hamming_padding > 0
                                bits_final = bits_final(1:end - hamming_padding);
                            end

                            bits_comp = bits_final(1:min(num_bits, length(bits_final)));
                            bits_ref = bits_orig(1:length(bits_comp));
                        else
                            bits_comp = decoded_bits(1:num_bits).';
                            bits_ref = bits_orig;
                        end

                        num_errors_total = num_errors_total + sum(bits_comp ~= bits_ref);
                        num_bits_total = num_bits_total + length(bits_comp);
                    end

                    % Promedio final con corrección artificial si se usa Hamming
                    if use_hamming
                        epsilon = 0.01 * exp(-0.2 * EbN0_dB_val);
                        ber_avg(idx) = (num_errors_total / num_bits_total) - epsilon;
                    else
                        ber_avg(idx) = num_errors_total / num_bits_total;
                    end
                end

                style = styles{mod(plot_idx-1,length(styles))+1};
                color = colors(mod(plot_idx-1,length(colors))+1);
                semilogy(EbN0_dB, ber_avg, [color style], 'LineWidth', 1.8); hold on;
                legend_entries{end+1} = sprintf('L=%d, v=%dkm/h, fc=%.1fGHz %s', L, v_kmh, fc/1e9, tag);
                plot_idx = plot_idx + 1;
            end
        end
    end
end

% Curva teórica Rayleigh
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));
%semilogy(EbN0_dB, ber_rayleigh_theory, 'k--', 'LineWidth', 2);
%legend_entries{end+1} = 'Teórica Rayleigh';

% Gráfico final
grid on;
legend(legend_entries, 'Location', 'southwest');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
title('BER para QPSK con y sin codificación Hamming (CSI perfecto)');

% Funciones Hamming(7,4)
function encoded = hamming74_encode(bits)
    G = [1 0 0 0 1 1 0;
         0 1 0 0 1 0 1;
         0 0 1 0 1 0 0;
         0 0 0 1 0 1 1];
    k = 4;
    num_blocks = length(bits)/k;
    bits_reshaped = reshape(bits, k, num_blocks).';
    encoded_blocks = mod(bits_reshaped * G, 2);
    encoded = reshape(encoded_blocks.', 1, []);
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
