function p = raised_cosine_fir_filter(T, t, alfa)
% Raised-Cosine FIR filter
    p = (sin(pi*t/T)./(pi*t/T)) .* (cos(alfa*pi*t/T)./(1-(2*alfa*t/T).^2));
end

