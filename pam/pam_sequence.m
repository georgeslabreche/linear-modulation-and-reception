function pam = pam_sequence()
% PAM sequence
    % PAM symbols
    a = [-3 -1 1 3]; % symbol alphabet/constellation.
    ind = randi(4, 100, 1); % random vector that includes integers between 1 and 4.
    pam = a(ind); % a sequence of random 4-PAM symbols.
end

