function psk_det = symbol_detection(snpsk, const_psk)
%symbol_detection 
%   n symbol detection, we compare the distance of a received noisy sample 
%   to all possible symbol values in the used alphabet.
%   The detected symbol (decision) is the one that minimizes the distance.

    % Nx20000-matrix, each line contains the same sn-vector.
    sn_block = repmat(snpsk, 1, length(const_psk)).';

    % Nx20000-matrix, where each column contains const.
    const_block = repmat(const_psk, 1, 20000);

    % Nx20000-matrix, whose every column contains the received symbol
    % distances to all possible symbol constellation points.
    distance = abs(sn_block - const_block); 

    % returns the minimum distance y and the corresponding
    % constellation index ind_1. Both vectors have the size of 1x20000.
    [~, ind] = min(distance); 

    % using vector ind, we can determine the detected symbol vector.
    psk_det = const_psk(ind);

end

