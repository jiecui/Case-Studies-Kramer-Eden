% %% [markdown]
% # Chapter 11 Analysis of Spike-Field Coherence during Navigation

% %% [markdown]
% ## Step 1 Visual inspection

% %% [markdown]
% ### Load the data

% %%
load('Ch11-spikes-LFP-1.mat') % load the multiscale dat

% %% [markdown]
% ### Plot the first trail

% %%
figure
plot(t, y(1, :))
hold on
plot(t, n(1, :), 'k')
hold off
ylim([-1, 1] * 1.2)
xlabel('Time [s]')
title('Figure 11.1 LFP and spikes in the first trial')

% %% [markdown]
% ### Spike-triggered average (STA)

% %%
[t_sta, STA] = spike_triggered_average(n, y, t);

% %% [markdown]
% #### Plot STA

% %%
figure
subplot(2, 1, 1)
plot(t_sta, STA([1, 6, 10, 16], :)')
xlim tight
ylim([-1, 1] * 0.4)
hold on
plot([0, 0], ylim, 'r')
legend(["trial 1", "trial 6", "trial 10", "trial 16"])
xlabel('Time [s]')

subplot(2, 1, 2)
plot(t_sta, STA', 'k')
xlim tight
ylim([-1, 1] * 0.25)
hold on
plot(t_sta, mean(STA), 'c', linewidth = 2)
plot([0, 0], ylim, 'r')
xlabel('Time [s]')

sgtitle('Figure 11.2 STA for initial trials (top) and all trials (bottom)')

% %% [markdown]
% ### Field-triggered average (FTA)

% %%
load('Ch11-spikes-LFP-1.mat') % load the multiscale dat

% 9-11 Hz
bp1 = [9, 11];
[FTA1, p] = field_triggered_average(n, y, t, bp1);

% 44-46 Hz
bp2 = [44, 46];
FTA2 = field_triggered_average(n, y, t, bp2);

% %% [markdown]
% #### Plot FTA

% %%
figure
subplot(1, 2, 1)
plot(p, mean(FTA1, 1))
xlim([-pi, pi])
ylim([0, .3])
xlabel('Phase [rad]')
ylabel('PTA')
title('9-11 Hz')

subplot(1, 2, 2)
plot(p, mean(FTA2, 1))
xlim([-pi, pi])
ylim([0, .3])
xlabel('Phase [rad]')
ylabel('PTA')
title('44-46 Hz')

sgtitle('Figure 11.3 FTA for LFP')

% %% [markdown]
% ## Step 2 Spike-field coherency

% %% [markdown]
% ### Add Chronux paths

% %%
addpath(genpath('~/Documents/Richard/ComputationalToolbox/neurophysiology_signals_analysis/chronux'))

% %% [markdown]
% ### Load the data

% %%
load('Ch11-spikes-LFP-1.mat') % load the multiscale dat

% %% [markdown]
% ### Estimate spke-field coherence using Chronux toolbox

% %%
N = size(y, 2);
dt = t(2) - t(1);
Fs = 1 / dt;
TW = 3;
ntapers = 2 * TW - 1;

% set the parameters of the MTM
params.Fs = Fs;
params.tapers = [TW, ntapers];
params.pad = -1;
params.trialave = 1;

% compute the MTM cohernence
[C, ~, ~, Syy, Snn, f] = coherencycpb(y.', n.', params);

% %% [markdown]
% #### Plot SFC

% %%
lambda = mean(sum(n, 2)) / (N * dt);

figure
subplot(1, 3, 1)
plot(f, pow2db(Snn))
hold on
plot(xlim, [1 1] * pow2db(lambda), 'r')
xlim([0, 100])
ylim([17, 23])
xlabel('Frequency [Hz]')
ylabel('Power density [dB/Hz]')
title('Snn')

subplot(1, 3, 2)
plot(f, pow2db(Syy))
xlim([0, 100])
ylim([-45, -15])
xlabel('Frequency [Hz]')
ylabel('Power density [dB/Hz]')
title('Syy')

subplot(1, 3, 3)
plot(f, C)
xlim([0, 100])
ylim([0, 0.5])
xlabel('Frequency [Hz]')
ylabel('Coherence')
title('Spike-field coh')

sgtitle('Figure 11.4 Coherence of LFP and spikes using MTM')

% %% [markdown]
% ### Rescale LFP

% %%
y_rs = .1 * y;
C_rs = coherencycpb(y_rs.', n.', params);

figure
plot(f, C_rs, f, C, '.')
legend({'Rescaled', 'Original'})
xlim([0, 100])
ylim([0 .5])
xlabel('Frequency [Hz]')
ylabel('Coherence')
title('Figure 11.5 SFC for scaled and original LFP')

% %% [markdown]
% ### Rescale spike trains by thinning

% %%
load('Ch11-spikes-LFP-1.mat') % load the multiscale dat
dt = t(2) - t(1);
N = size(y, 2);

thinning_factor = 0:.05:.75; %Choose a thinning factor.
num_tf = length(thinning_factor);
lambda = zeros(1, num_tf);
C_tf = zeros(num_tf, length(f));

for j = 1:num_tf
    tf_j = thinning_factor(j);
    n_j = n;

    for k = 1:size(n_j, 1) %For each trial,
        spikes = find(n_j(k, :) == 1); %...find the spikes.
        n_spikes = length(spikes); %...determine % of spikes.
        spikes = spikes(randperm(n_spikes)); %...permute spikes indices,
        n_remove = floor(tf_j * n_spikes); %...% spikes to remove,
        n_j(k, spikes(1:1 + n_remove)) = 0; %... remove the spikes.
    end % for

    lambda(j) = mean(sum(n_j, 2)) / (N * dt);
    C_tf(j, :) = coherencycpb(y.', n_j.', params);

end % for

% %%
figure
subplot(1, 2, 1)
plot(thinning_factor, lambda, '.-')
axis tight
axis square
xlabel('Thinning factor')
ylabel('Expected spike rate [Hz]')

subplot(1, 2, 2)
idx = ismember(thinning_factor, [0, .25, .5, .75]);
plot(f, C_tf(idx, :).')
legend(["0", "0.25", "0.5", "0.75"])
xlim([35, 55])
ylim([0, .5])
xlabel('Frequency [Hz]')
ylabel('Coherence')

sgtitle('Figure 11.6 SFC with original and thinned spike trains')

% %% [markdown]
% ## Step 3 Model of SFC using GLM

% %% [markdown]
% ### Estimate the model parameters

% %%
load('Ch11-spikes-LFP-1.mat') % load the multiscale dat

phi0 = -pi:.01:pi;

% 9-11 Hz
pb1 = [9, 11]; %Set the passband,
[pta1, phase1, phi1] = field_triggered_average(n, y, t, pb1);
[md1, ev1] = spike_field_glm(n, phi1, EvalPhase = phi0);

fprintf('9-11 Hz, beta_1 p-value: %f, beta_2 p-value: %f\n', ...
    md1.stats.p(2), md1.stats.p(3))

% 44-46 Hz
pb2 = [44, 46]; %Set the passband,
[pta2, phase2, phi2] = field_triggered_average(n, y, t, pb2);
[md2, ev2] = spike_field_glm(n, phi2, EvalPhase = phi0);

fprintf('44-46 Hz, beta_1 p-value: %f, beta_2 p-value: %f\n', ...
    md2.stats.p(2), md2.stats.p(3))

% %%
co = colororder();

figure
subplot(1, 2, 1)
y0 = ev1.y0;
dylo = ev1.dylo;
dyhi = ev1.dyhi;

plot(phase1, mean(pta1, 1), 'Color', co(1, :))
hold on
plot(phi0, y0, 'Color', co(2, :), 'LineWidth', 2)
plot(phi0, y0 + dylo, 'Color', co(2, :), 'LineWidth', 1)
plot(phi0, y0 - dyhi, 'Color', co(2, :), 'LineWidth', 1)
xlim tight
ylim([0 .25])
xlabel('Phase [rad]')
ylabel('Probability of a spike')
title('9-11 Hz')

subplot(1, 2, 2)
y0 = ev2.y0;
dylo = ev2.dylo;
dyhi = ev2.dyhi;

plot(phase2, mean(pta2, 1), 'Color', co(1, :))
hold on
plot(phi0, y0, 'Color', co(2, :), 'LineWidth', 2)
plot(phi0, y0 + dylo, 'Color', co(2, :), 'LineWidth', 1)
plot(phi0, y0 - dyhi, 'Color', co(2, :), 'LineWidth', 1)
xlim tight
ylim([0 .25])
xlabel('Phase [rad]')
ylabel('Probability of a spike')
title('44-46 Hz')

sgtitle('Figure 11.7 FTA and GLM fit')

% %% [markdown]
% ### Nested model

% %%
load('Ch11-spikes-LFP-1.mat') % load the multiscale dat

phi = ones(size(n));
X0 = phi(:);
Y = n(:);

[b0, dev0, stats0] = glmfit(X0, Y, 'poisson', 'constant', 'off');
dev2 = md2.dev;
pval = 1 - chi2cdf(dev0 - dev2, 2);
fprintf('p-value: %f\n', pval)

% %% [markdown]
% ### Thinning factor and GLM model

% %%
alpha = 0.05;
pd = makedist('normal');
Z = pd.icdf(1 - alpha / 2);

N = size(y, 2);
pb = [44, 46];

thinning_factor = 0:.05:.75; %Choose a thinning factor.
num_tf = length(thinning_factor);
fp = zeros(1, num_tf); % firing probability
beta_0 = zeros(1, num_tf); % model of overall firing probability
beta_1 = zeros(1, num_tf);
beta_2 = zeros(1, num_tf);

ci_0 = zeros(2, num_tf); % standard errors of beta
ci_1 = zeros(2, num_tf);
ci_2 = zeros(2, num_tf);

for j = 1:num_tf
    tf_j = thinning_factor(j);
    n_j = n;

    for k = 1:size(n_j, 1) %For each trial,
        spikes = find(n_j(k, :) == 1); %...find the spikes.
        n_spikes = length(spikes); %...determine % of spikes.
        spikes = spikes(randperm(n_spikes)); %...permute spikes indices,
        n_remove = floor(tf_j * n_spikes); %...% spikes to remove,
        n_j(k, spikes(1:1 + n_remove)) = 0; %... remove the spikes.
    end % for

    fp(j) = mean(sum(n_j, 2)) / N;

    % fit glm model
    [~, ~, phi_j] = field_triggered_average(n_j, y, t, pb);
    mod_j = spike_field_glm(n_j, phi_j);
    beta_0(j) = mod_j.b(1);
    beta_1(j) = mod_j.b(2);
    beta_2(j) = mod_j.b(3);

    ci_0(:, j) = mod_j.b(1) + [-1, 1] * Z * mod_j.stats.se(1);
    ci_1(:, j) = mod_j.b(2) + [-1, 1] * Z * mod_j.stats.se(2);
    ci_2(:, j) = mod_j.b(3) + [-1, 1] * Z * mod_j.stats.se(3);
end % for

% %%
co = colororder();

figure

subplot(1, 2, 1)
plot(thinning_factor, fp, 'Color', co(1, :), 'LineWidth', 2)
hold on
plot(thinning_factor, exp(beta_0), 'Color', co(2, :), 'LineWidth', 2)
plot(thinning_factor, exp(ci_0(1, :)), 'Color', co(2, :), 'LineStyle', ':')
plot(thinning_factor, exp(ci_0(2, :)), 'Color', co(2, :), 'LineStyle', ':')
axis tight square
xlabel('Thinning factor')
ylabel('Firing probability')

subplot(1, 2, 2)
plot(thinning_factor, exp(beta_1), 'Color', co(1, :), 'LineWidth', 2)
hold on
plot(thinning_factor, exp(ci_1(1, :)), 'Color', co(1, :), 'LineStyle', ':')
plot(thinning_factor, exp(ci_1(2, :)), 'Color', co(1, :), 'LineStyle', ':')
plot(thinning_factor, exp(beta_2), 'Color', co(2, :), 'LineWidth', 2)
plot(thinning_factor, exp(ci_2(1, :)), 'Color', co(2, :), 'LineStyle', ':')
plot(thinning_factor, exp(ci_2(2, :)), 'Color', co(2, :), 'LineStyle', ':')
xlim tight
ylim([.9 1.5])
plot(xlim, [1 1], 'g')
axis square
xlabel('Thinning factor')
ylabel('Exponentiated parameters')

sgtitle('Figure 11.8 GLM model and thinning factor: 44-46 Hz')


