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

%% 2. Modulation (SSB-SC AM)
%% Time Domain
% Hilbert Transfrom
mh1 = 2 * sin(2 * pi * f1 .* t);
mh2 = 2 * sin(2 * pi * f2 .* t);
mh3 = 2 * sin(2 * pi * f3 .* t);
mh = mh1 + mh2 + mh3;

% Modulate signals (Upper SSB AM signal)
u1 = m1 .* cos(2 * pi * fc .* t)-mh1.*sin(2*pi*fc.*t);
u2 = m2 .* cos(2 * pi * fc .* t)-mh2.*sin(2*pi*fc.*t);
u3 = m3 .* cos(2 * pi * fc .* t)-mh3.*sin(2*pi*fc.*t);
u = u1 + u2 + u3;
% Plot the figures
% Total period
figure('Name','modulated signal total');
plot(t,u);
% 1 period
figure('Name','modulated signal 1period');
plot(t(1:p),u(1:p));

%% Frequency Domain
% Perform the Fourier transform
U = fft(u) / length(u);
U_shift = fftshift(U);

% Plot the figures
% Total frequencies
figure('Name','modulated signal freq.domain');
plot(f, abs(U_shift));
% Partial frequencies
range2 = 1600;
index2 = (fs / 2 + 1) - range2 : (fs / 2 + 1) + range2;
f_part2 = f(index2);
U_shift_part = U_shift(index2);
figure('Name','partial freq domain modulated signal');
plot(f_part2, abs(U_shift_part));

%% 3. Demodulation (SSB-SC AM)
%% Time Domain
% deModulate signals
du = u .* cos(2 * pi * fc .* t);


%% Frequency Domain
% Perform the Fourier transform
DU = fft(du) / length(du);
DU_shift = fftshift(DU);
range5 = 2500;
index5 = (fs / 2 + 1) - range5 : (fs / 2 + 1) + range5;
f_part5 = f(index5);
DU_shift_part = DU_shift(index5);
% Plot the figures
% Total frequencies
figure();
plot(f_part5, abs(DU_shift_part));

%product signal through an ideal lowpass filter
range3=1000;
index3=(fs / 2 + 1) - range3:(fs / 2 + 1) + range3;
% BW=fc 범위의 저역통과 필터를 만들어 곱해준다.(사각펄스 신호)
LPF=1:fs;
LPF(fs/2+1-500:fs/2+1+500)=1;
LPF(1:fs/2+1-500)=0;
LPF(fs/2+1+500:fs)=0;
LPF_DU=DU_shift.* LPF;
figure('Name','Y(f) spectrum');
plot(f, abs(LPF_DU));

f_part3 = f(index3);
LPF_DU_shift_part = LPF_DU(index3);
figure();
plot(f_part3, abs(LPF_DU_shift_part));
%역 퓨리에 변환 과정=>y(t)
yshift=ifftshift(LPF_DU);
yt=ifft(yshift)*length(yshift);
figure('Name','y(t) total period');
plot(t,yt);
figure('Name','y(t) 1period');
plot(t(1:p),yt(1:p));

%% Modulation VSB Filter
% product carrier signal
uc1 = m1 .* cos(2 * pi * fc .* t);
uc2 = m2 .* cos(2 * pi * fc .* t);
uc3 = m3 .* cos(2 * pi * fc .* t);
uc = uc1 + uc2 + uc3;
% FT
UC = fft(uc) / length(uc);
UC_shift = fftshift(UC);

% define H(f): 사다리꼴 모양의 필터
fa=300;
W=2000;
H=zeros(1, fs);
H(fs/2+1-(fc+W):fs/2+1-fc-fa)=1;
H(fs/2+1+(fc+fa):fs/2+1+fc+W)=1;
for n = fs/2+1-(fc+fa):fs/2+1-(fc-fa)
        H(n) = 1 - (n-(fs/2+1-fc-fa)) / (2*fa);
end

for n = fs/2+1+(fc-fa):fs/2+1+(fc+fa)
        H(n) = (n-(fs/2+1+fc-fa)) / (2*fa);
end

for n = fs/2+1-(fc+500):fs/2+1-(fc+W)
        H(n) = (n-(fs/2+1-fc-500)) / (500-W);
end

for n = fs/2+1+(fc+W):fs/2+1+(fc+500)
        H(n) = 1-(n-(fs/2+1+fc+W)) / (500-W);
end

figure('Name','H(f) sidebandfilter');
plot(fs/2+1-4000:fs/2+1+4000,H(fs/2+1-4000:fs/2+1+4000));

UC_shift_part=UC_shift(index2);
figure();
plot(f_part2,abs(UC_shift_part));

SBF_U=UC_shift .*H;
figure('Name','amp spectrum of U(f) after pass SBF');
SBF_U_part=SBF_U(index2);
plot(f_part2, abs(SBF_U_part));

su_shift=ifftshift(SBF_U);
su=ifft(su_shift)*length(su_shift);

%% Demodulation
v=su .* cos(2*pi*fc.*t);

figure('Name','v(t) of 1 period');
plot(t(1:p),v(1:p));

V=fft(v)/length(v);
V_shift=fftshift(V);

figure();
plot(f,abs(V_shift));

%LPF 필터를 만든다.
LPF2=zeros(1,fs);
LPF2(fs/2+1-fc:fs/2+1+fc)=1;
for n=fs/2+1-fc-100:fs/2+1-fc
    LPF2(n)=(n-(fs/2+1-fc-100))/(100);
end
for n=fs/2+1+fc:fs/2+1+fc+100
    LPF2(n)=1-(n-(fs/2+1+fc))/(100);
end

figure();
plot(f,abs(LPF2));

LPF_V=V_shift .* LPF2;
LPF_V_part=LPF_V(index1);
figure('Name','V(f) after pass LPF');
plot(f,abs(LPF_V));
figure('Name','V(f) after pass LPF');
plot(f_part1,abs(LPF_V_part));

y2_shift=ifftshift(LPF_V);
y2=ifft(y2_shift)*length(y2_shift);

figure('Name','y(t)');
plot(t(1:p),y2(1:p));

figure('Name','result');
subplot(4,1,1);
plot(t(1:p),m(1:p));
title('message signal of 1period');

subplot(4,1,2);
plot(t(1:p),su(1:p));
title('modulated signal of 1period');

subplot(4,1,3);
plot(t(1:p),v(1:p));
title('demodulated signal of 1period before pass LPF');

subplot(4,1,4);
plot(t(1:p),y2(1:p));
title('demodulated signal of 1period after pass LPF');

figure('Name','result of frequency domain');
subplot(4,1,1);
plot(f_part1, abs(M_shift_part));
title('message signal');

subplot(4,1,2);
plot(f, abs(SBF_U));
title('modulated signal of 1period');

subplot(4,1,3);
plot(f, abs(V_shift));
title('demodulated signal before pass LPF');

subplot(4,1,4);
plot(f_part1, abs(LPF_V_part));
title('demodulated signal after pass LPF');