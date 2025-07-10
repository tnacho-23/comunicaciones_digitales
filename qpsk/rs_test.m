clc; clear all; close all;

% -------------------------------
% Parámetros generales
% -------------------------------
use_rs_coding = true;  % <--- Activar/desactivar codificación RS
num_bits = 1e4;        % Cantidad base de bits (se ajustará si se usa RS)
k = 2;                 % bits por símbolo QPSK
EbN0_dB = -2:2:30;
num_runs = 21;
mapping = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 0; 1 1];
snr_plot_vals = [-2, 0, 10, 30];
saved_constellations = struct();
L_vals = [10];
v_kmh_vals = [50];          
fc_vals = [700e6, 3.5e9];        
colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y', 'c'];
styles = {'-', '--', '-.', ':'};

% -------------------------------
% Reed-Solomon Parameters
% -------------------------------
if use_rs_coding
    m = 8;
    n_rs = 2^m - 1;
    k_rs = 239;
    bits_per_rsword = k_rs * m;
    rsEncoder = comm.RSEncoder('CodewordLength', n_rs, 'MessageLength', k_rs);
    rsDecoder = comm.RSDecoder('CodewordLength', n_rs, 'MessageLength', k_rs);

    min_bits_rs = bits_per_rsword;
    if num_bits < min_bits_rs
        error('num_bits (%d) es muy pequeño para aplicar RS. Usa al menos %d bits.', num_bits, min_bits_rs);
    end

    num_rswords = floor(num_bits / bits_per_rsword);
    num_bits = num_rswords * bits_per_rsword;
    fprintf('Ajustando num_bits a %d para RS (múltiplo de %d bits)\n', num_bits, bits_per_rsword);
end

figure;
legend_entries = {};
plot_idx = 1;

% -------------------------------
% Bucle principal con Doppler
% -------------------------------
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
            bits_orig = randi([0 1], 1, num_bits);

            if use_rs_coding
                bits_rs_in = reshape(bits_orig, 8, []).';
                sym_rs = bi2de(bits_rs_in, 'left-msb');
                sym_enc = step(rsEncoder, sym_rs);
                bits_enc = de2bi(sym_enc, 8, 'left-msb');
                bits = reshape(bits_enc.', 1, []);
            else
                bits = bits_orig;
            end

            num_symbols = length(bits) / k;

            bit_pairs = reshape(bits, 2, []).';
            idx_sym = bit_pairs(:,1)*2 + bit_pairs(:,2) + 1;
            symbols = mapping(idx_sym).';

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

            decoded_bits = zeros(length(bits),1);
            for n = 1:num_symbols
                d = abs(y_eq(n) - mapping).^2;
                [~, idx_min] = min(d);
                decoded_bits(2*n-1:2*n) = bit_map(idx_min,:).';
            end

            if use_rs_coding
                bits_dec = reshape(decoded_bits, 8, []).';
                sym_rec = bi2de(bits_dec, 'left-msb');
                sym_corr = step(rsDecoder, sym_rec);
                bits_corr = de2bi(sym_corr, 8, 'left-msb');
                bits_final = reshape(bits_corr.', 1, []);
                ber_total(idx) = ber_total(idx) + sum(bits_final ~= bits_orig(1:length(bits_final)));
            else
                ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);
            end

            if run == 1 && iL==1 && iv==1 && ifc==1 && ismember(EbN0_dB(idx), snr_plot_vals)
                if EbN0_dB(idx) < 0, s = sprintf('snr_m%d', abs(EbN0_dB(idx)));
                else, s = sprintf('snr_%d', EbN0_dB(idx)); end
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

% -------------------------------
% Caso sin Doppler
% -------------------------------
ber_total = zeros(1, length(EbN0_dB));
for idx = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(idx)/10);
    N0 = 1/(2*k*EbN0);
    for run = 1:num_runs
        bits_orig = randi([0 1], 1, num_bits);

        if use_rs_coding
            bits_rs_in = reshape(bits_orig, 8, []).';
            sym_rs = bi2de(bits_rs_in, 'left-msb');
            sym_enc = step(rsEncoder, sym_rs);
            bits_enc = de2bi(sym_enc, 8, 'left-msb');
            bits = reshape(bits_enc.', 1, []);
        else
            bits = bits_orig;
        end

        num_symbols = length(bits) / k;

        bit_pairs = reshape(bits, 2, []).';
        idx_sym = bit_pairs(:,1)*2 + bit_pairs(:,2) + 1;
        symbols = mapping(idx_sym).';

        H = (randn(1, num_symbols) + 1j*randn(1, num_symbols)) / sqrt(2);
        noise = sqrt(N0)*(randn(1,num_symbols) + 1j*randn(1,num_symbols));
        y = symbols .* H + noise;
        y_eq = y ./ H;

        decoded_bits = zeros(length(bits),1);
        for n = 1:num_symbols
            d = abs(y_eq(n) - mapping).^2;
            [~, idx_min] = min(d);
            decoded_bits(2*n-1:2*n) = bit_map(idx_min,:).';
        end

        if use_rs_coding
            bits_dec = reshape(decoded_bits, 8, []).';
            sym_rec = bi2de(bits_dec, 'left-msb');
            sym_corr = step(rsDecoder, sym_rec);
            bits_corr = de2bi(sym_corr, 8, 'left-msb');
            bits_final = reshape(bits_corr.', 1, []);
            ber_total(idx) = ber_total(idx) + sum(bits_final ~= bits_orig(1:length(bits_final)));
        else
            ber_total(idx) = ber_total(idx) + sum(decoded_bits.' ~= bits);
        end
    end
end

ber_avg = ber_total / (num_runs * num_bits);
semilogy(EbN0_dB, ber_avg, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);
legend_entries{end+1} = sprintf('L=inf');

% -------------------------------
% Curva teórica Rayleigh
% -------------------------------
EbN0_lin = 10.^(EbN0_dB/10);
ber_rayleigh_theory = 0.5*(1 - sqrt(EbN0_lin./(EbN0_lin+1)));
semilogy(EbN0_dB, ber_rayleigh_theory, 'k--', 'LineWidth', 2);
legend_entries{end+1} = 'Teórica Rayleigh';

grid on;
legend(legend_entries, 'Location', 'southwest');
xlabel('E_b/N_0 [dB]');
ylabel('BER');
title('BER para QPSK en canal Multipath + AWGN');

% -------------------------------
% Mostrar constelaciones
% -------------------------------
figure;
snrs = snr_plot_vals;
titles = {'1. Original', '2. Solo AWGN', '3. Multipath + AWGN', '4. Ecualizado'};

for i = 1:length(snrs)
    if snrs(i) < 0, s = sprintf('snr_m%d', abs(snrs(i)));
    else, s = sprintf('snr_%d', snrs(i)); end

    if isfield(saved_constellations, s)
        data = saved_constellations.(s);
        subplot(4,4,(i-1)*4+1); plot(real(data.tx), imag(data.tx), 'bo'); axis equal; grid on; title(titles{1});
        ylabel(sprintf('SNR = %d dB', snrs(i)));
        subplot(4,4,(i-1)*4+2); plot(real(data.rx_awgn), imag(data.rx_awgn), 'ro'); axis equal; grid on; title(titles{2});
        subplot(4,4,(i-1)*4+3); plot(real(data.rx_multipath), imag(data.rx_multipath), 'mo'); axis equal; grid on; title(titles{3});
        subplot(4,4,(i-1)*4+4); plot(real(data.rx_eq), imag(data.rx_eq), 'go'); axis equal; grid on; title(titles{4});
    else
        fprintf('Advertencia: no se guardaron constelaciones para SNR = %d dB\n', snrs(i));
    end
end
sgtitle('Constelaciones para SNR = -2, 0, 10, 30 dB (L=5, v=30km/h, fc=700MHz)');
