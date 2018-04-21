addpath('./pam/');
addpath('./qsk/');

clear all 
close all
clc

%--------------------%
% 1.1 4-PAM Sequence %
%--------------------%
pam = pam_sequence();

% plot
figure;
stem(pam);

%--------------------%
% 1.2 Nyquist Pulse  %
%--------------------%
Fs = 24000; % Sampling frequency 24000 Hz
T = 1/8000; % Symbol time interval [s].
t = -5*T:1/Fs:5*T; % Time vector (sampling intervals)
t = t+1e-10; % Otherwise, the denominator would be zero at t=0
% (or manually set p (t=0)=1 )
alfa = 0.5; % Roll-off factor

p = raised_cosine_fir_filter(T, t, alfa);

% plot
figure;
plot(t,p);
hold on;
stem(t,p);
xlabel('Time [s]');
ylabel('Amplitude');
hold off;

%--------------------------------------------%
% 1.3 Linear modulation with a Nyquist pulse %
%--------------------------------------------%
r = Fs*T; % Oversampling factor
xn = pulse_shaping_filtering(pam, p, r);

% plot
figure;
plot(xn(1:200)); % A piece of the signal

%---------------------%
% 1.4 The Eye-Diagram %
%---------------------%
figure
hold on;
for i = 16:6:291
    plot(xn(i:i+6));
end
hold off;
grid on;

var_zn = 0.01; % Set noise variance
zn = sqrt(var_zn)*randn(size(xn)); % Generate noise with defined variance
yn = xn + zn; % Add noise to the signal

% Plot the eye-diagram for yn as above
figure
hold on;
for i = 16:6:291
    plot(yn(i:i+6));
end
hold off;
grid on;


%---------------------------------%
% 1.4 The Special Characteristics %
%---------------------------------%
Nfft = 2048; % Size of the FFT
P = fftshift(fft(p,Nfft)); % Spectrum (FFT) of the pulse-shape
X = fftshift(fft(xn,Nfft)); % Spectrum (FFT) of the PAM waveform
Y = fftshift(fft(yn,Nfft)); % Spectrum (FFT) of the noisy PAM waveform
f = -Fs/2:Fs/Nfft:Fs/2-Fs/Nfft; % Frequency axis for plotting the spectra

% Figure: Pulse Spectrum.
figure;
subplot(211);
plot(f,abs(P));
grid;
title('Pulse Spectrum');
xlabel('Frequency [Hz]');
subplot(212);
plot(f,20*log10(abs(P)));
grid;
title('Pulse Spectrum in dB');
xlabel('Frequency [Hz]');

% Figure: Signal Spectrum in dB.
figure;
subplot(211);
plot(f,abs(X));
grid;
title('Signal Spectrum');
xlabel('Frequency [Hz]');
subplot(212);
plot(f,20*log10(abs(X)));
grid;
title('Signal Spectrum in dB');
xlabel('Frequency [Hz]');

% Figure: Noisy Signal Spectrum.
figure;
subplot(211);
plot(f,abs(Y));
grid;
title('Noisy Signal Spectrum');
xlabel('Frequency [Hz]');
subplot(212);
plot(f,20*log10(abs(Y)));
grid;
title('Noisy Signal Spectrum in dB');
xlabel('Frequency [Hz]');

