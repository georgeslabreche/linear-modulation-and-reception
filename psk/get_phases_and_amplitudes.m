function result = get_phases_and_amplitudes(const_psk, psk_sequence)
%get_phases_and_amplitudes
    % Returns an object with the following properties
    %   - Number of phases.
    %   - List of all phases in radians (unique values).
    %   - List of all phases in degrees (unique values).
    %   - Number of amplitudes.
    %   - List of different amplitude values (unique values).
    
    % Set values we are interested in.
    result.phases.rad = unique(angle(const_psk));
    result.phases.deg = radtodeg(result.phases.rad);
    result.phase_count = length(result.phases.rad);
    result.amplitudes = unique(abs(psk_sequence));
    result.amplitude_count = length(result.amplitudes);
end

