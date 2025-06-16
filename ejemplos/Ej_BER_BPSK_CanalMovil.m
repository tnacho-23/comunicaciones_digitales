clc
clear all
tic
M = 1e5;    %número de bits/símbolos
b = round( rand(1,M) );
d = 2*b-1;                  %No-retorno cero BPSK, modulación antípoda
SNRdB = 0:1:5;

%señal modulada interferente
b_i = round( rand(1,M) );
d_i = b_i-0.5;   

%---------------------------------------------------------------------
%parametros del canal multipath considerando la física del sistema
 fc = 2.4e9; %portadora central
 lambda = 3e8/fc; %longitud de onda
 v = (50)/3.6; %velocidad del móvil UE
 fmax = v/lambda; %(fmax = v/lambda)
%Dt = 10e-8;
%t = (0:0.1:M-1)*1;
t = linspace(0,1,M);
N =10; %número de paths del ganal
an = ones(1,N)*sqrt(1/N); %potencia normalizada de los paths
thetan = 2*pi*rand(1,N); %fase aleatoria del path
fDn = fmax * cos( 2*pi*rand(1,N) ); %Doppler aleatorio
H = zeros(1,M);
 
 for n = 1:N
    H = H + an(n)*exp( j*( thetan(n) - 2*pi*fDn(n)*t ) );
 end
 
%-------------------------------------------------------------- 
%canal Rayleigh multipath multiplicativo (N=infinito)%
%H = ( randn( 1, M ) + j*randn( 1, M ) )*sqrt(1/2);

H_i = ( randn( 1, M ) + j*randn( 1, M ) )*sqrt(1/2); 
%Canal Rayleigh
%considerando señal interferente co-canal


%%
Pe = [];
for n = 1:length(SNRdB)
    constante=rand(1);
    Pn = 10^( -SNRdB(n)/10 );
    n = ( randn( 1, M ) + j*randn( 1, M ) )*sqrt(Pn/2);
    y = (d.*H)+n+(d_i.*H_i); %Interferencia
    %y = (d.*H)+n;
    yh = y./H; %ecualización, se asume estimación perfecta del canal
    z = sign( real(yh) ); %detección de la señal antípoda
    Pe = [Pe length( find( z-d ) )/M]; %errrores del sistema
end

figure(1)
subplot(2,1,1)
plot(t, 10*log10(abs(H).^2));

subplot(2,1,2)
plot(t, angle(H));

figure(2)
SNR = 10.^(SNRdB/10);
PeAWGN = 0.5*erfc( sqrt(SNR) );
PeT = 0.5*(1-sqrt( SNR./(1+SNR) ));
semilogy( SNRdB, PeAWGN, 'k' )
hold on
semilogy( SNRdB, Pe, 'b' )
semilogy( SNRdB, PeT, '--r' )
hold off
grid on
axis( [0 SNRdB(end) 1e-6 1] )
legend('BER en un canal AWGN', 'BER en un canal movil: Simulacion', 'BER en un canal movil: Teoria')
toc

figure(3)
subplot(3,1,1)
plot(real(d),imag(d),'x');
title('BPSK PLOT');
xlabel('REAL(d)');
ylabel('IMG(d)');
legend('Símbolos Transmitidos','PLOT AT TX' );

subplot(3,1,2)
plot(real(d+n),imag(d+n),'x');
title('BPSK PLOT');
xlabel('REAL(d)');
ylabel('IMG(d)');
legend('Símbolos Transmitidos Con Ruido','PLOT AT TX' );

subplot(3,1,2)
plot(real(d+n),imag(d+n),'x');
title('BPSK PLOT');
xlabel('REAL(d)');
ylabel('IMG(d)');
legend('Símbolos Transmitidos Con Ruido','PLOT AT TX' );

subplot(3,1,3)
plot(real(y),imag(y),'x');
title('BPSK PLOT');
xlabel('REAL(d)');
ylabel('IMG(d)');
legend('Símbolos Recibidos con canal Rayleigh y AWGN sin Ecualizar','PLOT AT TX' );

figure(4)

subplot(4,1,1)
plot(real(y.*H),imag(y.*H),'x');
title('BPSK PLOT');
xlabel('REAL(d)');
ylabel('IMG(d)');
legend('Símbolos Recibidos con canal Rayleigh  sin Ecualizar','PLOT AT TX' );

subplot(4,1,2)
plot(real(y),imag(y),'x');
title('BPSK PLOT');
xlabel('REAL(d)');
ylabel('IMG(d)');
legend('Símbolos Recibidos con canal Rayleigh + AWGN sin Ecualizar','PLOT AT TX' );


subplot(4,1,3)
plot(real(yh),imag(yh),'x');
title('BPSK PLOT');
xlabel('REAL(d)');
ylabel('IMG(d)');
legend('Símbolos Recibidos Ecualizado','PLOT AT TX' );