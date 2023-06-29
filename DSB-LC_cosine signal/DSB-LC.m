close all
clear
clc

%% simulation spec
fs=8192;
t=0:1/fs:1-1/fs;
f1=100;
f2=200;
f3=300;
fc=1000;

%% 1.Message signal
%% Time Domain
% define message signal
m1 = 2 * cos(2 * pi * f1 .* t);
m2 = 2 * cos(2 * pi * f2 .* t);
m3 = 2 * cos(2 * pi * f3 .* t);
m = m1 + m2 + m3;

% Plot the figures
% Total period
figure();
plot(t, m1);
hold on
plot(t, m2);
plot(t, m3);
plot(t, m);
% 1 period
p = ceil(fs / f1); % fs / f1 = 1 period
figure('Name','m1,m2,m3,m graph(Time Domain)');
plot(t(1 : p), m1(1 : p));
hold on
plot(t(1 : p), m2(1 : p));
plot(t(1 : p), m3(1 : p));
plot(t(1 : p), m(1 : p));
%% Frequency Domain
% Perform the Fourier transform
M = fft(m) / length(m);
M_shift = fftshift(M);

% Plot the figures
% Total frequencies
f = -fs / 2 : fs / 2 - 1;
figure();
plot(f, abs(M_shift));
% Partial frequencies
range1 = 400;
index1 = (fs / 2 + 1) - range1 : (fs / 2 + 1) + range1;
f_part1 = f(index1);
M_shift_part = M_shift(index1);
figure('Name','M graph');
plot(f_part1, abs(M_shift_part));

%% 2. Modulation (DSB-LC AM)
%% Time Domain
%mn = m/abs(m)  abs(m)=6
mn=m/max(m);
a=0.5; 
%modulation index
u = (a*mn+1).*cos(2*pi*fc.*t);
figure();
plot(t,u);
figure('Name','u time domain grph');
plot(t(1 : p), u(1 : p));
%% Frequency Domain
% Perform the Fourier transform
U = fft(u) / length(u);
U_shift = fftshift(U);

% Plot the figures
% Total frequencies
figure();
plot(f, abs(U_shift));
% Partial frequencies
range2 = 1500;
index2 = (fs / 2 + 1) - range2 : (fs / 2 + 1) + range2;
f_part2 = f(index2);
U_shift_part = U_shift(index2);
figure('Name','U freq domain grph');
plot(f_part2, abs(U_shift_part));

%% Demodulate (DSB-LC)
%#1 using envelope detector
%rectify only positive signal
up=u>=0;     %find what index has positive value, and that index has true=1 other has false=0
uup=u.*up;  %then we get only positive value by producting u(t) with up
figure('Name','uup');
plot(t(1:p),uup(1:p));
figure('Name','enveope detector time domian');
plot(t(1:p),uup(1:p));
DU=fft(uup) / length(uup);
DU_shift=fftshift(DU);
%demodulate using envelope detector


range3=1000;
index3=(fs / 2 + 1) - range3:(fs / 2 + 1) + range3;
% #2 BW=500 Hz 범위의 저역통과 필터를 만들어 곱해준다.(사각펄스 신호)
LPF=1:fs;
LPF(fs/2+1-500:fs/2+1+500)=1;
LPF(1:fs/2+1-500)=0;
LPF(fs/2+1+500:fs)=0;
%passing through to LPF
LPF_DU=DU_shift.* LPF;
figure('Name','D(f) before passing LPF');
plot(f,abs(DU_shift));
figure('Name','D(f) passing LPF');
plot(f, abs(LPF_DU));

f_part3 = f(index3);
LPF_DU_shift_part = LPF_DU(index3);
figure('Name','product LPF filter');
plot(f_part3, abs(LPF_DU_shift_part));

%역퓨리에 변환 과정=>d(t)
LPF_DU_ishift=ifftshift(LPF_DU);
d=ifft(LPF_DU_ishift)*length(LPF_DU_ishift);
figure('Name','1period envelope detector d(t) after LPF filter');
plot(t(1:p),d(1:p));
% d(t) 로부터 m(t)신호 추출
y=(d-abs(LPF_DU(fs/2+1)));  % substract DC component
Y=fft(y)/length(y);
Y_shift=fftshift(Y);
Y_shift_part=Y_shift(index1);
%Y(f) 출력
figure('Name','Y(f) is');
plot(f,abs(Y_shift));
%Y(f) 의 일부분만 출력
figure();
plot(f_part1,abs(Y_shift_part))
range4=3000;
index4=(fs / 2 + 1) - range4:(fs / 2 + 1) + range4;
f_part4=f(index4);
DU_shift_part=DU_shift(index4);

%%Plot all result of time domain signal
figure('Name','result of time domain signal');
subplot(4,1,1);
plot(t(1:p),m(1:p));
title('message signal of 1period');

subplot(4,1,2);
plot(t(1:p),u(1:p));
title('modulated signal of 1period');

subplot(4,1,3);
plot(t(1:p),uup(1:p));
title('output envelope detector signal of 1period upper signal');

subplot(4,1,4);
plot(t(1:p),y(1:p));
title('demodulated signal y(t) of 1period');

%% Plot all result of frequency domain signal
figure('Name','result of frequency domain');
subplot(4,1,1);
plot(f_part1, abs(M_shift_part));
title('message signal');

subplot(4,1,2);
plot(f_part2, abs(U_shift_part));
title('modulated signal of 1period');

subplot(4,1,3);
plot(f_part4, abs(DU_shift_part));
title('demodulated signal before pass LPF');

subplot(4,1,4);
plot(f_part1, abs(Y_shift_part));
title('demodulated signal after pass LPF (DC component delete and scailing)');

