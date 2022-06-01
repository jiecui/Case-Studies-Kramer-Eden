function [model, eval] = spike_field_glm(n, phi, options)
    % GLM model of spike-field interactions/association

    arguments
        n (:, :) double {mustBeNumeric, mustBeReal}
        phi (:, :) double {mustBeNumeric, mustBeReal} % instantaneous phase
    end % positional

    arguments
        options.EvalPhase = -pi:0.01:pi;
    end % optional

    eval_pha = options.EvalPhase;

    phi = phi(:); %Convert phase matrix to vector.
    X = [cos(phi) sin(phi)]; %Create a matrix of predictors.
    Y = n(:); %Create a vector of responses.
    [b, dev, stats] = glmfit(X, Y, 'poisson'); %Fit the GLM.

    phi0 = transpose(eval_pha); %Define new phase interval,
    X0 = [cos(phi0) sin(phi0)]; %...and predictors,
    [y0, dylo, dyhi] = glmval(b, X0, 'log', stats); %...evaluate the model.

    % ouput
    model.b = b;
    model.dev = dev;
    model.stats = stats;

    eval.y0 = y0;
    eval.dylo = dylo;
    eval.dyhi = dyhi;
end % function

% [EOF]
