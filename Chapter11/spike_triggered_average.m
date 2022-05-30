function [t, sta] = spike_triggered_average(spikes, lfp, time, options)

    arguments
        spikes (:, :) {mustBeNumeric}
        lfp (:, :) {mustBeNumeric}
        time(1, :) double
    end % positional

    arguments
        options.WindowLength (1, 1) double = 100 % in samples
    end % optional

    win = options.WindowLength;
    n = spikes;
    y = lfp;
    dt = time(2) - time(1);
    t = linspace(-win * dt, win * dt, 2 * win + 1);

    K = size(n, 1);
    N = size(y, 2);
    sta = zeros(K, 2 * win + 1);

    for k = 1:K
        spks = find(n(k, :) == 1);

        for i = 1:length(spks)

            if spks(i) > win && spks(i) < N - win
                sta(k, :) = sta(k, :) + y(k, spks(i) - win:spks(i) + win) / length(spks);
            end % if

        end % for

    end % for

end % function

% [EOF]
