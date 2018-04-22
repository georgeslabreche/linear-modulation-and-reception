addpath('./psk/');

clear all 
close all
clc

% Show plots or just export them directly as an image file (or both!).
show_plots = 'on';
export_plots = false;

% Create plot image export directory if it doesn't exist.
export_dir = 'plots/psk/';
fn = fullfile(export_dir);
if ~exist(fn, 'dir') && export_plots
   mkdir(export_dir);
end

%----------------------------------------------------%
% 2.1 The symbol constellation and a symbol sequence %
%----------------------------------------------------%
% QPSK % 
% QPSK is a special case of M-PSK (for M = 4).
% In M-PSK, the constellation points are located
% at regular angular intervals on the unit circle.
const_qpsk = exp(1j*[pi/4 3*pi/4 5*pi/4 7*pi/4]).'; % QPSK alphabet.
qpsk = const_qpsk(randi(4,20000,1)); % QPSK symbol sequence.

% Get the phase and amplitudes.
pcaa_qpsk = get_phases_and_amplitudes(const_qpsk, qpsk)
pcaa_qpsk_rads = pcaa_qpsk.phases.rad
pcaa_qpsk_degs = pcaa_qpsk.phases.deg

% Plot the QPSK (4-PSK) constellation.
f = figure('Name', 'QPSK Constellation');
set(f, 'Visible', show_plots);
plot(qpsk,'o');
xlabel('Re');
ylabel('Im');
title('QPSK Constellation');
axis([-2 2 -2 2]); % scale the axis of the figure.
if export_plots == true
    print(strcat(export_dir, 'qpsk-constellation.png'), '-dpng');
end

% 8-PSK %
% (Q18) Also generate an 8-PSK symbol stream and plot the constellation.
% 8-PSK alphabet.
const_8psk = exp(1j*[0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4]).'; 
psk8 = const_8psk(randi(8,20000,1)); % 8-PSK symbol sequence.

% Get the phase and amplitudes.
pcaa_8psk = get_phases_and_amplitudes(const_8psk, psk8)
pcaa_8psk_rads = pcaa_8psk.phases.rad
pcaa_8psk_degs = pcaa_8psk.phases.deg

% Plot the 8-PSK constellation.
f = figure('Name', '8-PSK Constellation');
set(f, 'Visible', show_plots);
plot(psk8,'o');
xlabel('Re');
ylabel('Im');
title('8-PSK Constellation');
axis([-2 2 -2 2]); % scale the axis of the figure.
if export_plots == true
    print(strcat(export_dir, '8-psk-constellation.png'), '-dpng');
end

% 64-QAM % 
aqam = [-7 -5 -3 -1 1 3 5 7];
A = repmat(aqam,8,1);
B = flipud(A') ;
const_qam = A+1j*B; % 8x8-matrix with constellation points.
const_qam = const_qam(:); % column-vector with 64-QAM alphabet.
qam = const_qam(randi(64,20000,1)); % 64-QAM symbol sequence.

% Get the phase and amplitudes.
pcaa_qam = get_phases_and_amplitudes(const_qam, qam)
pcaa_qam_rads = pcaa_qam.phases.rad
pcaa_qam_degs = pcaa_qam.phases.deg
pcaa_qam_amplitudes = pcaa_qam.amplitudes

% Plot 64-QAM  constellation.
f = figure('Name', '64-QAM Constellation');
set(f, 'Visible', show_plots);
plot(qam,'o');
axis([-10 10 -10 10]); % Scale the axis of the figure.
xlabel('Re');
ylabel('Im');
title('64-QAM Constellation');
if export_plots == true
    print(strcat(export_dir, '64-qam-constellation.png'), '-dpng');
end

%--------------------------------%
% 2.2 Adding Noise to the Signal %
%--------------------------------%
n = randn(size(qam))+1j*randn(size(qam)); % randn generates Gaussian noise.
svqpsk = std(qpsk)^2; % QPSK signal variance (power).
sv8psk = std(psk8)^2; % 8-PSK signal variance (power).
svqam = std(qam)^2; % QAM signal variance (power).
nv = std(n)^2; % noise variance (power).

SNR = 20; % TargetSNR=20dB

% add noise to the QPSK signal.
pqpsk = std(qpsk)/(std(n)*10^(SNR/20)); % proper constant p.
snqpsk = qpsk + n * pqpsk; % add noise to signal.

% plot noisy QPSK constellation.
f = figure('Name', 'Noisy QPSK Constellation');
set(f, 'Visible', show_plots);
plot(snqpsk,'o');
xlabel('Re');
ylabel('Im');
title('Noisy QPSK Constellation');
axis([-2 2 -2 2]); % scale the axis of the figure.
if export_plots == true
    print(strcat(export_dir, 'noisy-qpsk-constellation.png'), '-dpng');
end

% add noise to the 8-PSK signal.
p8psk = std(psk8)/(std(n)*10^(SNR/20)); % proper constant p. 
sn8psk = psk8 + n * p8psk; 

% plot noisy 8-PSK constellation.
f = figure('Name', 'Noisy 8-PSK Constellation');
set(f, 'Visible', show_plots);
plot(sn8psk,'o');
xlabel('Re');
ylabel('Im');
title('Noisy 8-PSK Constellation');
axis([-2 2 -2 2]); % scale the axis of the figure.
if export_plots == true
    print(strcat(export_dir, 'noisy-8-psk-constellation.png'), '-dpng');
end

% add noise to 64-QAM signal.
pqam = std(qam)/(std(n)*10^(SNR/20)); % proper constant p.
snqam = qam + n * pqam; 

% plot noisy 64-QAM constellation.
f = figure('Name', 'Noisy 64-QAM Constellation');
set(f, 'Visible', show_plots);
plot(snqam,'o');
axis([-10 10 -10 10]); % Scale the axis of the figure.
xlabel('Re');
ylabel('Im');
title('Noisy 64-QAM Constellation');
if export_plots == true
    print(strcat(export_dir, 'noisy-64-qam-constellation.png'), '-dpng');
end

%---------------------------------------------------%
% 2.3 Signal Detection (symbol by symbol detection) %
%---------------------------------------------------%
%------%
% QPSK %
%------%
% 4x20000-matrix, each line contains the same sn-vector.
sn_block = repmat(snqpsk,1,4).';

% 4x20000-matrix, where each column contains const.
const_block = repmat(const_qpsk,1,20000);

% 4x20000-matrix, whose every column contains the received symbol
% distances to all possible symbol constellation points.
distance = abs(sn_block-const_block); 

% returns the minimum distance y and the corresponding
% constellation index ind_1. Both vectors have the size of 1x20000.
[~,ind_1] = min(distance); 

% using vector ind_1, we can determine the detected symbol vector.
qpsk_det = const_qpsk(ind_1);

%-----%
% QAM %
%-----%
% 64x20000-matrix, each line contains the same sn-vector.
sn_block = repmat(snqam,1,64).';

% 64x20000-matrix, where each column contains the const.
const_block = repmat(const_qam,1,20000);

% 64x20000-matrix, whose every column contains the received symbol 
% distances to all possible symbol constellation points.
distance = abs(sn_block-const_block);

% returns the minimum distance y and the corresponding
% constellation index ind_2. Both vectors have the size of 1x20000.
[y,ind_2] = min(distance); 

% using vector ind_2, we can determine the detected symbol vector.
qam_det = const_qam(ind_2); 


%---------------------------------%
% 2.4 The symbol-error rate (SER) %
%---------------------------------%
% minimum distance d for our QPSK constellation.
d = sqrt(2);

% sigma is the deviation of noise (real or imaginary part).
sigma = std(real(n * pqpsk));

Q = 0.5*erfc(d/(sqrt(2)*2*sigma));

% theoretical symbol error rate.
ser_theo = 2*Q - Q^2

% Comparison returns 0 (false) or 1 (true)
nr_of_errors = sum(qpsk~=qpsk_det);

% simulated symbol error rate.
ser_simu = nr_of_errors/20000

d = 2; % minimum distance d for our QAM alphabet

% sigma is the deviation of noise real or imaginary part
sigma = std(real(n * pqam));
Q = 0.5*erfc(d/(sqrt(2)*2*sigma));

nr_of_errors = sum(qam~=qam_det);
ser_simu = nr_of_errors/20000 % simulated symbol error rate
ser_theo = 3.5*Q - 3.0625*Q^2 % theoretical symbol error rate