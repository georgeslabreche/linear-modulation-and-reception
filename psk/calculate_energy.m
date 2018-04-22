function energy = calculate_energy(psk)
%CALCULATE_ENERGY
    energy = sum(real(psk).^2 + imag(psk).^2);
end

