function [fta, phase, phi] = field_triggered_average(n, y, t, pass_band, options)

    arguments
        n (:, :) double {mustBeNumeric, mustBeReal}
        y (:, :) double {mustBeNumeric, mustBeReal}
        t (1, :) double {mustBeNumeric, mustBeReal} % in seconds
        pass_band (1, 2) double {mustBeNumeric, mustBeReal} % in Hz
    end % positional

    arguments
        options.FilterOrder (1, 1) double = 100
    end % optional

    % parameters
    ord = options.FilterOrder;
    dt = t(2) - t(1);
    fNQ = 1 / dt / 2;

    K = size(n, 1);
    N = size(y, 2);
    Wn = pass_band / fNQ;
    b = fir1(ord, Wn);

    % estimate fta
    phase = linspace(-pi, pi, N);
    fta = zeros(K, N);
    phi = zeros(K, N); %Create variable to hold phase.

    for k = 1:K
        Vlo = filtfilt(b, 1, y(k, :));
        phi_k = angle(hilbert(Vlo)); % instantaneous phase of signal
        [~, indices] = sort(phi_k);
        fta(k, :) = n(k, indices);
        phi(k, :) = phi_k;
    end % for

end % function

% [EOF]
