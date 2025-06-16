clear
%% H(t) = sum_(n=1)^N an exp( j( thetan - 2*pi*fDn*t ) )
%Modelo de canal de Banda Angosta
clc
clear all
fc = 2.4e9; %portadora
lambda = 3e8/fc; %longitud de onda de la portadora
v = (250)/3.6; %velocidad del UE movil
fmax = v/lambda; %(fmax = v/lambda) %frecuencia Doppler máxima
M = 1e5;    %cantidad de muestras
t = linspace(0,0.25,M);


%-----------------------------------------------------------------------
N =10000;  %paths
an = ones(1,N)*sqrt(1/N); %energía normalizada de los paths
thetan = 2*pi*rand(1,N);    %fase aleatoria
fDn = fmax * cos( 2*pi*rand(1,N) );     %Frecuencia Doppler dependiente del angulo de arrivo.

H = zeros(1,M);
for n = 1:N
    H = H + an(n)*exp( j*( thetan(n) - 2*pi*fDn(n)*t ) );
end
%-------------------------------------------------------------------------
%Por el teorema del límite central, el canal puede modelarse como dos
%procesos gaussianos (real e imaginario). Se comporta como un ruido
%multiplicativo. 
%H = (randn(1,M) +j*randn(1,M))*sqrt(1/2);
% Y = X*H + N
%--------------------------------------------------------------------
%Plot de magnitud y fase de la función de transferencia del canal
figure(1)
subplot(2,1,1)
plot(t, 10*log10(abs(H).^2));

subplot(2,1,2)
plot(t, angle(H));

%-----------------------------------------------------------------
%PDF PARTE REAL E IMAGINARIA
figure(2)
z = linspace(-5,5,100);
Dz = abs( z(1)-z(2) );
pdfAI = hist(imag(H),z)/(M*Dz); %PDF IMAGINARIO
plot(z,pdfAI)
hold on

pdfAR = hist(real(H),z)/(M*Dz); %PDF REAL
plot(z,pdfAR)
sigma=sqrt(0.5);
pdfAT = (1/(sqrt(2*pi)*sigma))*exp( -z.^2/(2*sigma^2) );
hold on
plot(z,pdfAT)
hold off
%------------------------------------------------------------
%FASE
figure(3)
theta = linspace(-1.5*pi, 1.5*pi,100);
Dtheta = abs( theta(1)-theta(2) );
pdfT = hist(angle(H),theta)/(M*Dtheta);
plot(theta,pdfT)
