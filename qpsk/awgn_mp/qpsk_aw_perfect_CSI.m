clc; clear all; close all;

% Parámetros fijos
num_bits = 1e4;
k = 2;
num_symbols = num_bits / k;
EbN0_dB = -2:2:30;
num_runs = 10;  
mapping = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 0; 1 1];

% Variables a recorrer
L_vals = [5, 40];
v_kmh_vals = [30, 120];
fc_vals = [700e6, 3.5e9];

% Colores y estilos para distinguir curvas
colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y', 'c'];  % Al menos 8 colores
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

            ber_total = zeros(1, length(EbN0_dB));

            for idx = 1:length(EbN0_dB)
                EbN0 = 10^(EbN0_dB(idx)/10);
                N0 = 1/(2*k*EbN0);

                for run = 1:num_runs
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

            % Plot
            style = styles{mod(plot_idx-1,length(styles))+1};
            color = colors(mod(plot_idx-1,length(colors))+1);
            semilogy(EbN0_dB, ber_avg, [color style], 'LineWidth', 1.8); hold on;

            % Etiqueta
            legend_entries{end+1} = sprintf('L=%d, v=%dkm/h, fc=%.1fGHz', L, v_kmh, fc/1e9);
            plot_idx = plot_idx + 1;
        end
    end
end

% Teórica
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));
semilogy(EbN0_dB, ber_rayleigh_theory, 'k--', 'LineWidth', 2);

legend_entries{end+1} = 'Teórica Rayleigh';

grid on;
legend(legend_entries, 'Location', 'southwest');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
title('BER para QPSK en canal Multipath + AWGN');
