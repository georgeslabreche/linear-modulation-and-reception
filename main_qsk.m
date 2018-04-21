%addpath('./qsk/');

clear all 
close all
clc

%----------------------------------------------------%
% 2.1 The symbol constellation and a symbol sequence %
%----------------------------------------------------%
const_qpsk = exp(1j*[pi/4 3*pi/4 5*pi/4 7*pi/4]).'; % QPSK alphabet.
qpsk = const_qpsk(randi(4,20000,1)); % QPSK symbol sequence.

% plot the constellation.
figure
plot(qpsk,'o');
xlabel('Re');
ylabel('Im');
axis([-2 2 -2 2]); % scale the axis of the figure.

aqam = [-7 -5 -3 -1 1 3 5 7];
A = repmat(aqam,8,1);
B = flipud(A') ;
const_qam = A+1j*B; % 8x8-matrix with constellation points.
const_qam = const_qam(:); % column-vector with 64-QAM alphabet.
qam = const_qam(randi(64,20000,1)); % 64-QAM symbol sequence.

% plot the constellation.
figure
plot(qam,'o');
axis([-8 8 -8 8]); % Scale the axis of the figure.
xlabel('Re');
ylabel('Im');

%--------------------------------%
% 2.2 Adding Noise to the Signal %
%--------------------------------%
n = randn(size(qam))+1j*randn(size(qam)); % randn generates Gaussian noise.
svqpsk = std(qpsk)^2; % QPSK signal variance (power).
svqam = std(qam)^2; % QAM signal variance (power).
nv = std(n)^2; % noise variance (power).

SNR = 20; % TargetSNR=20dB

% add noise to signal.
p1 = std(qpsk)/(std(n)*10^(SNR/20)); % proper constant p.
snqpsk = qpsk+n*p1; % add noise to signal.

% plot noisy constellation.
figure
plot(snqpsk,'o'); 

% add noise to signal.
p2 = std(qam)/(std(n)*10^(SNR/20)); % proper constant p.
snqam = qam+n*p2; 

% plot noisy constellation.
figure
plot(snqam,'o');

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
sigma = std(real(n*p1));

Q = 0.5*erfc(d/(sqrt(2)*2*sigma));

% theoretical symbol error rate.
ser_theo = 2*Q - Q^2

% Comparison returns 0 (false) or 1 (true)
nr_of_errors = sum(qpsk~=qpsk_det);

% simulated symbol error rate.
ser_simu = nr_of_errors/20000

d = 2; % minimum distance d for our QAM alphabet

% sigma is the deviation of noise real or imaginary part
sigma = std(real(n*p2));
Q = 0.5*erfc(d/(sqrt(2)*2*sigma));

nr_of_errors = sum(qam~=qam_det);
ser_simu = nr_of_errors/20000 % simulated symbol error rate
ser_theo = 3.5*Q - 3.0625*Q^2 % theoretical symbol error rate