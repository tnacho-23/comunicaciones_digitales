clc; clear all; close all;

% Parámetros fijos
num_bits = 1e5;
k = 4;  % 4 bits por símbolo para 16-QAM
num_symbols = num_bits / k;
EbN0_dB = -2:2:30;
num_runs = 10;

% Mapeo 16-QAM (Gray coding)
real_vals = [-3 -1 3 1];
imag_vals = [-3 -1 3 1];
[re, im] = meshgrid(real_vals, imag_vals);
mapping = (re(:) + 1j*im(:)) / sqrt(10); % Normalizado a potencia unitaria

% Mapa de bits Gray correspondiente
bit_map = [
    0 0 0 0;
    0 0 0 1;
    0 0 1 1;
    0 0 1 0;
    0 1 1 0;
    0 1 1 1;
    0 1 0 1;
    0 1 0 0;
    1 1 0 0;
    1 1 0 1;
    1 1 1 1;
    1 1 1 0;
    1 0 1 0;
    1 0 1 1;
    1 0 0 1;
    1 0 0 0
];

% Valores de SNR para los que se graficarán constelaciones
snr_plot_vals = [-2, 0, 10, 30];
saved_constellations = struct();

% Variables a recorrer
L_vals = [5, 40];
v_kmh_vals = [30, 120];          
fc_vals = [700e6, 3.5e9];        

% Colores y estilos
colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y', 'c'];
styles = {'-', '--', '-.', ':'};

figure;
legend_entries = {};
plot_idx = 1;

% ----------------------------
% CASOS CON DOPPLER
% ----------------------------
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
                EbN0_dB_val = EbN0_dB(idx);
                EbN0 = 10^(EbN0_dB_val/10);
                N0 = 1/(2*k*EbN0);

                for run = 1:num_runs
                    bits = randi([0 1], 1, num_bits);
                    bit_groups = reshape(bits, k, []).';
                    
                    % Buscar el índice del símbolo en la tabla
                    indices = zeros(num_symbols, 1);
                    for i = 1:num_symbols
                        for j = 1:size(bit_map, 1)
                            if isequal(bit_groups(i, :), bit_map(j, :))
                                indices(i) = j;
                                break;
                            end
                        end
                    end
                    symbols = mapping(indices).';

                    % Canal y ruido con Doppler
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
                        decoded_bits((n-1)*k+1:n*k) = bit_map(idx_min,:).';
                    end
                    ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);

                    should_save = (run == 1) && (iL == 1) && (iv == 1) && (ifc == 1);
                    if should_save && ismember(EbN0_dB_val, snr_plot_vals)
                        if EbN0_dB_val < 0
                            s = sprintf('snr_m%d', abs(EbN0_dB_val));
                        else
                            s = sprintf('snr_%d', EbN0_dB_val);
                        end
                        saved_constellations.(s).tx = symbols;
                        saved_constellations.(s).rx_awgn = symbols + noise;
                        saved_constellations.(s).rx_multipath = y;
                        saved_constellations.(s).rx_eq = y_eq;
                        fprintf('Guardando constelaciones para %s\n', s);
                    end
                end
            end

            ber_avg = ber_total / (num_runs * num_bits);

            style = styles{mod(plot_idx-1,length(styles))+1};
            color = colors(mod(plot_idx-1,length(colors))+1);
            semilogy(EbN0_dB, ber_avg, [color style], 'LineWidth', 1.8); hold on;

            legend_entries{end+1} = sprintf('L=%d, v=%dkm/h, fc=%.1fGHz', L, v_kmh, fc/1e9);
            plot_idx = plot_idx + 1;
        end
    end
end

% ----------------------------
% CASO SIN DOPPLER
L_inf = 40;
color = [0.5 0.5 0.5];  
style = '--';

for iL = 1:length(L_inf)
    L = L_inf(iL);
    ber_total = zeros(1, length(EbN0_dB));

    for idx = 1:length(EbN0_dB)
        EbN0_dB_val = EbN0_dB(idx);
        EbN0 = 10^(EbN0_dB_val/10);
        N0 = 1/(2*k*EbN0);

        for run = 1:num_runs
            bits = randi([0 1], 1, num_bits);
            bit_groups = reshape(bits, k, []).';

            indices = zeros(num_symbols, 1);
            for i = 1:num_symbols
                for j = 1:size(bit_map, 1)
                    if isequal(bit_groups(i, :), bit_map(j, :))
                        indices(i) = j;
                        break;
                    end
                end
            end
            symbols = mapping(indices).';

            H = (randn(1, num_symbols) + 1j*randn(1, num_symbols)) / sqrt(2);
            noise = sqrt(N0)*(randn(1,num_symbols) + 1j*randn(1,num_symbols));
            y = symbols .* H + noise;
            y_eq = y ./ H;

            decoded_bits = zeros(num_bits,1);
            for n = 1:num_symbols
                distances = abs(y_eq(n) - mapping).^2;
                [~, idx_min] = min(distances);
                decoded_bits((n-1)*k+1:n*k) = bit_map(idx_min,:).';
            end
            ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);
        end
    end

    ber_avg = ber_total / (num_runs * num_bits);

    semilogy(EbN0_dB, ber_avg, 'Color', color, 'LineStyle', style, 'LineWidth', 2); hold on;
    legend_entries{end+1} = sprintf('L=inf');
end

% Teórica Rayleigh para 16-QAM (aproximada)
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 3/8*(1 - sqrt(EbN0_lin./(EbN0_lin+5)));  % Aproximación
semilogy(EbN0_dB, ber_rayleigh_theory, 'k--', 'LineWidth', 2);
legend_entries{end+1} = 'Teórica Rayleigh';

grid on;
legend(legend_entries, 'Location', 'southwest');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
title('BER para 16-QAM en canal Multipath + AWGN');

% Constelaciones
figure;
snrs = snr_plot_vals;
titles = {'1. Original', '2. Solo AWGN', '3. Multipath + AWGN', '4. Ecualizado'};

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

sgtitle('Constelaciones para SNR = -2, 0, 10, 30 dB (L=5, v=30km/h, fc=700MHz)');
