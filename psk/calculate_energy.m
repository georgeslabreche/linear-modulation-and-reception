function energy = calculate_energy(psk)
%calculate_energy
    energy = sum(real(psk).^2 + imag(psk).^2);
end

