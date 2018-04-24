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

energy_qpsk = calculate_energy(const_qpsk)

% Get the phase and amplitudes.
pcaa_qpsk = get_phases_and_amplitudes(const_qpsk, qpsk)
pcaa_qpsk_rads = pcaa_qpsk.phases.rad
pcaa_qpsk_degs = pcaa_qpsk.phases.deg

% Plot the QPSK (4-PSK) constellation.
f = figure('Name', 'QPSK Constellation');
set(f, 'Visible', show_plots);
plot(qpsk,'o');
xlabel('Real');
ylabel('Imaginary');
title('QPSK Constellation');
axis([-2 2 -2 2]); % scale the axis of the figure.
if export_plots == true
    print(strcat(export_dir, 'qpsk-constellation.png'), '-dpng');
end

% 8-PSK %
const_8psk = exp(1j*[0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4]).'; 
psk8 = const_8psk(randi(8,20000,1)); % 8-PSK symbol sequence.

% Calculate the energy.
energy_8psk = calculate_energy(const_8psk)

% Get the phase and amplitudes.
pcaa_8psk = get_phases_and_amplitudes(const_8psk, psk8)
pcaa_8psk_rads = pcaa_8psk.phases.rad
pcaa_8psk_degs = pcaa_8psk.phases.deg

% Plot the 8-PSK constellation.
f = figure('Name', '8-PSK Constellation');
set(f, 'Visible', show_plots);
plot(psk8,'o');
xlabel('Real');
ylabel('Imaginary');
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

% Calculate the energy.
energy_qam = calculate_energy(const_qam)

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
xlabel('Real');
ylabel('Imaginary');
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

for SNR = [20, 15, 10, 5]

    % add noise to the QPSK signal.
    pqpsk = std(qpsk)/(std(n)*10^(SNR/20)); % proper constant p.
    snqpsk = qpsk + n * pqpsk; % add noise to signal.

    % plot noisy QPSK constellation.
    plot_title = strcat('Noisy QPSK Constellation (', num2str(SNR), ' dB)');
    f = figure('Name', plot_title);
    set(f, 'Visible', show_plots);
    plot(snqpsk,'o');
    xlabel('Real');
    ylabel('Imaginary');
    %title(plot_title);
    axis([-2 2 -2 2]); % scale the axis of the figure.
    if export_plots == true
        print(strcat(export_dir, num2str(SNR), '-db-noise-qpsk-constellation.png'), '-dpng');
    end

    % add noise to the 8-PSK signal.
    p8psk = std(psk8)/(std(n)*10^(SNR/20)); % proper constant p. 
    sn8psk = psk8 + n * p8psk; 

    % plot noisy 8-PSK constellation.
    plot_title = strcat('Noisy 8-PSK Constellation (', num2str(SNR), ' dB)');
    f = figure('Name', plot_title);
    set(f, 'Visible', show_plots);
    plot(sn8psk,'o');
    xlabel('Real');
    ylabel('Imaginary');
    %title(plot_title);
    axis([-2 2 -2 2]); % scale the axis of the figure.
    if export_plots == true
        print(strcat(export_dir, num2str(SNR), '-db-noise-8-psk-constellation.png'), '-dpng');
    end

    % add noise to 64-QAM signal.
    pqam = std(qam)/(std(n)*10^(SNR/20)); % proper constant p.
    snqam = qam + n * pqam; 

    % plot noisy 64-QAM constellation.
    plot_title = strcat('Noisy 64-QAM Constellation (', num2str(SNR), ' dB)');
    f = figure('Name', plot_title);
    set(f, 'Visible', show_plots);
    plot(snqam,'o');
    axis([-10 10 -10 10]); % Scale the axis of the figure.
    xlabel('Real');
    ylabel('Imaginary');
    %title(plot_title);
    if export_plots == true
        print(strcat(export_dir, num2str(SNR), '-db-noise-64-qam-constellation.png'), '-dpng');
    end

end

%---------------------------------------------------%
% 2.3 Signal Detection (symbol by symbol detection) %
%---------------------------------------------------%
ser_theo_qpsk = [];
ser_simu_qpsk = [];
ser_simu_qam = [];
ser_theo_qam = [];
SNR_range = 0:25;

for SNR = SNR_range
    % Adding noise to QPSK and 64-QAM signals.
    % add noise to the QPSK signal.
    pqpsk = std(qpsk)/(std(n)*10^(SNR/20)); % proper constant p.
    snqpsk = qpsk + n * pqpsk; % add noise to signal.
    
    % add noise to 64-QAM signal.
    pqam = std(qam)/(std(n)*10^(SNR/20)); % proper constant p.
    snqam = qam + n * pqam; 

    % QPSK symbol detection.
    % using vector ind_1, we can determine the detected symbol vector.
    qpsk_det = symbol_detection(snqpsk, const_qpsk);

    % QAM symbol detection.
    % using vector ind_2, we can determine the detected symbol vector.
    qam_det = symbol_detection(snqam, const_qam); 

    %---------------------------------%
    % 2.4 The symbol-error rate (SER) %
    %---------------------------------%
    %------%
    % QPSK %
    %------%
    % minimum distance d for our QPSK constellation.
    d = sqrt(2);

    % sigma is the deviation of noise (real or imaginary part).
    sigma = std(real(n * pqpsk));

    Q = 0.5*erfc(d/(sqrt(2)*2*sigma));

    % theoretical symbol error rate.
    ser_theo_qpsk = [ser_theo_qpsk; 2*Q - Q^2];

    % simulated symbol error rate.
    % comparison returns 0 (false) or 1 (true)
    nr_of_errors = sum(qpsk~=qpsk_det);
    ser_simu_qpsk = [ser_simu_qpsk; nr_of_errors/20000];

    %-----%
    % QAM %
    %-----%
    d = 2; % minimum distance d for our QAM alphabet

    % sigma is the deviation of noise real or imaginary part
    sigma = std(real(n * pqam));
    Q = 0.5*erfc(d/(sqrt(2)*2*sigma));
    
    % theoretical symbol error rate
    ser_theo_qam = [ser_theo_qam; 3.5*Q - 3.0625*Q^2];

    % simulated symbol error rate
    nr_of_errors = sum(qam~=qam_det);
    ser_simu_qam = [ser_simu_qam; nr_of_errors/20000];
    
end

% Plot QPSK SER %
ser_theo_qpsk
ser_simu_qpsk
ser_diff_qpsk = (abs(ser_theo_qpsk - ser_simu_qpsk) ./ ((ser_theo_qpsk + ser_simu_qpsk) / 2)) * 100

plot_title = 'Theoretical and simulated SER for QPSK with SNR from 0 to 10.';
f = figure('Name', plot_title);
set(f, 'Visible', show_plots);
plot(SNR_range(1:11), ser_theo_qpsk(1:11));
hold on;
plot(SNR_range(1:11), ser_simu_qpsk(1:11));
xlabel('SNR');
ylabel('SER');
%title(plot_title);
legend('Theoretical', 'Simulated');
hold off;
if export_plots == true
    print(strcat(export_dir,'qpsk-ser-noise-0-to-10.png'), '-dpng');
end

plot_title = 'Theoretical and simulated SER for QPSK with SNR from 15 to 25.';
f = figure('Name', plot_title);
set(f, 'Visible', show_plots);
plot(SNR_range(16:26), ser_theo_qpsk(16:26));
hold on;
plot(SNR_range(16:26), ser_simu_qpsk(16:26));
xlabel('SNR');
ylabel('SER');
%title(plot_title);
legend('Theoretical', 'Simulated');
hold off;
if export_plots == true
    print(strcat(export_dir,'qpsk-ser-noise-15-to-25.png'), '-dpng');
end

% Plot QAM SER %
ser_theo_qam
ser_simu_qam
ser_diff_qam = ((ser_theo_qam - ser_simu_qam) ./ ((ser_theo_qam + ser_simu_qam) / 2)) * 100

plot_title = 'Theoretical and simulated SER for 64-QAM with SNR from 0 to 10.';
f = figure('Name', plot_title);
set(f, 'Visible', show_plots);
plot(SNR_range(1:11), ser_theo_qam(1:11));
hold on;
plot(SNR_range(1:11), ser_simu_qam(1:11));
xlabel('SNR');
ylabel('SER');
%title(plot_title);
legend('Theoretical', 'Simulated');
hold off;
if export_plots == true
    print(strcat(export_dir, '64-qam-ser-noise-0-to-10.png'), '-dpng');
end

plot_title = 'Theoretical and simulated SER for 64-QAM with SNR from 15 to 25.';
f = figure('Name', plot_title);
set(f, 'Visible', show_plots);
plot(SNR_range(16:26), ser_theo_qam(16:26));
hold on;
plot(SNR_range(16:26), ser_simu_qam(16:26));
xlabel('SNR');
ylabel('SER');
%title(plot_title);
legend('Theoretical', 'Simulated');
hold off;
if export_plots == true
    print(strcat(export_dir,'64-qam-ser-noise-15-to-25.png'), '-dpng');
end
