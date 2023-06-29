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
figure();
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
figure();
plot(f_part1, abs(M_shift_part));

%% 2. Modulation (DSB-SC AM)
%% Time Domain
% Modulate signals
u1 = m1 .* cos(2 * pi * fc .* t);
u2 = m2 .* cos(2 * pi * fc .* t);
u3 = m3 .* cos(2 * pi * fc .* t);
u = u1 + u2 + u3;
% Plot the figures
% Total period
figure('Name','total period of u(t)');
plot(t,u);
% 1 period
figure('Name','1 period of u(t)');
plot(t(1:p),u(1:p));
%% Frequency Domain
% Perform the Fourier transform
U = fft(u) / length(u);
U_shift = fftshift(U);

% Plot the figures
% Total frequencies
figure('Name','amp spectrum graph of U(f)');
plot(f, abs(U_shift));
% Partial frequencies
range2 = 1500;
index2 = (fs / 2 + 1) - range2 : (fs / 2 + 1) + range2;
f_part2 = f(index2);
U_shift_part = U_shift(index2);
figure('Name','amp spectrum graph of U(f)');
plot(f_part2, abs(U_shift_part));

%% 3. Demodulation (DSB-SC AM)
%% Time Domain
% deModulate signals
du1 = u1 .* cos(2 * pi * fc .* t);
du2 = u2 .* cos(2 * pi * fc .* t);
du3 = u3 .* cos(2 * pi * fc .* t);
du = du1 + du2 + du3;
% Plot the figures
% Total period
figure('Name','total period of r(t)');
plot(t,du);
%1period
figure('Name','1period of r(t)');
plot(t(1:p),du(1:p));

%% Frequency Domain
% Perform the Fourier transform
DU = fft(du) / length(du);
DU_shift = fftshift(DU);
% Plot the figures
% Total frequencies
figure('Name','amp spectrum of Y(f) before pass LPF');
plot(f, abs(DU_shift));
% partial frequency
range5=2500;
index5=(fs / 2 + 1) - range5:(fs / 2 + 1) + range5;
f_part5 = f(index5);
DU_shift_part = DU_shift(index5);
figure('Name','amp spectrum of Y(f) before pass LPF');
plot(f_part5, abs(DU_shift_part));
%product signal through an ideal lowpass filter
range3=1000;
index3=(fs / 2 + 1) - range3:(fs / 2 + 1) + range3;
% BW=fc 범위의 저역통과 필터를 만들어 곱해준다.(사각펄스 신호)
LPF=zeros(1,fs);
LPF(fs/2+1-fc:fs/2+1+fc)=1;

%product LPF to demodulated signal
LPF_DU=DU_shift.* LPF;
figure('Name','amp spectrum of Y(f) after pass LPF');
plot(f, abs(LPF_DU));
f_part3 = f(index3);
%f_part4
range4=3000;
index4=(fs / 2 + 1) - range4:(fs / 2 + 1) + range4;
f_part4=f(index4);
DU_shift_part=DU_shift(index4);
LPF_DU_shift_part = LPF_DU(index1);
figure('Name','amp spectrum of Y(f) after pass LPF');
plot(f_part1, abs(LPF_DU_shift_part));

%역 퓨리에 변환 과정=>y(t)
y_shift=ifftshift(LPF_DU);
y=ifft(y_shift)*length(y_shift);
figure('Name','y(t)');
plot(t(1:p),y(1:p));

%% Plot all result signal of time domain
figure('Name','result');
subplot(4,1,1);
plot(t(1:p),m(1:p));
title('message signal of 1period');

subplot(4,1,2);
plot(t(1:p),u(1:p));
title('modulated signal of 1period');

subplot(4,1,3);
plot(t(1:p),du(1:p));
title('demodulated signal of 1period before pass LPF');

subplot(4,1,4);
plot(t(1:p),y(1:p));
title('demodulated signal of 1period after pass LPF');

%% Plot all result of frequency domain signal
figure('Name','result of frequency domain');
subplot(4,1,1);
plot(f_part1, abs(M_shift_part));
title('message signal');

subplot(4,1,2);
plot(f_part2, abs(U_shift_part));
title('modulated signal ');

subplot(4,1,3);
plot(f_part4, abs(DU_shift_part));
title('demodulated signal before pass LPF');

subplot(4,1,4);
plot(f_part1, abs(LPF_DU_shift_part));
title('demodulated signal after pass LPF');

