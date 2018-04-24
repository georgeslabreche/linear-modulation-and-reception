function diff = percentage_difference(a, b)
%percentage_difference
%   Calculates the percentage difference between two values.
    diff = (abs(a - b) ./ ((a + b) / 2)) * 100;
end

