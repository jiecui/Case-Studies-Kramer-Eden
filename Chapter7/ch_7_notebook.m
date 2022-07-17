% %% [markdown]
% # Chapter 7 Investigation of Cross-Frequency Coupling in a Local Field Potential

% %% [markdown]
% ## Step 1 Visual inspection

% %% [markdown]
% ### Load the data

% %%
load('Ch7-LFP-1.mat')

% %% [markdown]
% ### Plot LFP

% %%
figure
plot(t, LFP)
axis([4 5 -2 2])
hold on
plot(ones(2, 1) * [4.35, 4.5], ylim' * ones(1, 2), 'r--')
xlabel('Time [s]')
ylabel('Voltage [mV]')
title('Figure 7.1 Example trace of LFP')

% %% [markdown]
% ## Step 2 Estimate power spectrum

% %% [markdown]
% ### 2.1 Using classic FFT

% %%
dt=t(2)-t(1);
Fs=1/dt;
fNQ = Fs / 2;
T=t(end);
N=length(LFP);

% %% [markdown]
% #### estimate FFT using Hanning window

% %%
x = hann(N) .* transpose(LFP);
xf = fft(x - mean(x)); % zero-mean process
Sxx = 2 * dt^2 / T * (xf .* conj(xf));
Sxx = Sxx(1:N / 2 + 1);

% %% [markdown]
% #### plot the power spectrum

% %%
df = 1 / max(T);
faxis = (0:df:fNQ);
figure
plot(faxis, 10 * log10(Sxx))
xlim([0 200])
ylim([-80 0])
grid on
grid minor
xlabel('Frequency (Hz)')
ylabel('Power [mV^2/Hz]')
title('Figure 7.2 Spectrum of LFP data')

% %% [markdown]
% ### 2.2 Using multitaper method to estimate spectrum

% %% [markdown]
% #### add directories

% %% [markdown]
% * on cortex

% %%
addpath(genpath('~/Documents/Richard/Projects/rochester/mayo_clinic/nonuniform_signal'))

% %% [markdown]
% #### estimate spectrum

% %%
TW = 40;
f_qf = 0:.1:fNQ / 2;
nus = NUSpectrum(x - mean(x), t);
J = nus.mtnufft('QuerryFrequencies', f_qf, 'TimeHalfbandwidth', TW);
J = J / 3.5/10;

% %% [markdown]
% #### plot

% %%
pxx = mean(J .* conj(J), 2);
figure
h_fft = plot(faxis, pow2db(Sxx));
hold on
h_mt = plot(f_qf, pow2db(pxx), 'LineWidth', 3);
xlim([0 200])
ylim([-80 0])
grid on
grid minor
legend([h_fft, h_mt], ["FFT", "MTFFT"])
xlabel('Frequency (Hz)')
ylabel('Power [mV^2/Hz]')
title('Figure 7.2a MT spectrum of LFP data')

% %% [markdown]
% #### F-test of periodicity

% %%
Fval = nus.nuftest(QuerryFrequencies = f_qf, TimeHalfbandwidth = TW);
Nr = 2 * fNQ / df;
K = 2 * TW - 1;
sig = finv(1 - 1 / Nr, 2, 2 * K - 2);

figure
plot(f_qf, Fval)
hold on
plot(xlim, ones(1, 2) * sig, 'r')
xlim([0, 200])
xlabel('Frequency (Hz)')
ylabel('F-statistic')
title('Figure 7.2b F-test of LFP data')

% %% [markdown]
% ## Step 3 Cross-frequency coupling

% %% [markdown]
% ### Filtering the data into high- and low-frequency bands

% %% [markdown]
% #### for low-frequency interval

% %%
Wn = [5, 7] / fNQ; % set the passband
n = 100; % filter order
b = fir1(n, Wn);
Vlo = filtfilt(b, 1, LFP); % apply the filter

% %% [markdown]
% #### for high-frequency interval

% %%
Wn = [80, 120] / fNQ;
b = fir1(n, Wn);
Vhi = filtfilt(b, 1, LFP);

% %% [markdown]
% #### plot Figure 7.3

% %%
D = timetable(LFP', Vlo', Vhi', 'RowTimes', seconds(t'), ...
'VariableNames', ["LFP", "Low", "High"]); % make the timetable

% stacked plot
figure
stackedplot(D)
xlim(seconds([4 5]))

% %% [markdown]
% ### plot Figure 7.4 

% %%
phi = angle(hilbert(Vlo)); % compute phase of low frequency signal
amp = abs(hilbert(Vhi)); % compute amplitude of high frequency signal

% %% [markdown]
% #### plot Figure 7.4a

% %%
figure
ax_lo = subplot(2, 1, 1);
plot(t, Vlo)
title('Low frequency signal')

ax_hi = subplot(2, 1, 2);
plot(t, phi);
hold on
plot(xlim' * ones(1, 2), ones(2, 1) * [-1, 1] * pi, 'c')
xlabel('Time [s]')
title('Phase')

linkaxes([ax_lo, ax_hi], 'x')
xlim(ax_hi, [4 5])
sgtitle('Figure 7.4a Low frequency signal and phase')


% %% [markdown]
% #### plot Figure 7.4b

% %%
figure
plot(t, Vhi)
hold on
plot(t, amp, 'LineWidth', 2)
xlim([4 5])
xlabel('Time [s]')
title('Figure 7.4b High frequency signal and amplitude')


