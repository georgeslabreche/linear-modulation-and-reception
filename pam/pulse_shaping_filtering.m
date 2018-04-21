function xn = pulse_shaping_filtering(pam, p, r)
% Pulse shaping filtering
    N = length(pam); % Number of symbols
    pams = zeros(size(1:r*N));
    pams(1:r:r*N) = pam; % Up-sampled sequence
    xn = filter(p,1,pams); % Pulse shaping filtering
end

