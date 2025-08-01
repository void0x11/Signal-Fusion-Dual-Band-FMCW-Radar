% R2 FMCW Radar Model – Fixed Rx Length and Chirp Count Issues
clc; clear; close all;

%% === 1. Constants and Radar Specs ===
c = 3e8;
fc = 5.8e9;
lambda = c / fc;
Rmax = 100;                 % Maximum range
Ts = (2*Rmax / c);
bw = 150e6;                 % نكبر BW بنفس النسبة عشان الميل μ يظل ثابت
fs = 2 * bw;                % لازم fs يتوافق مع BW الجديد
sweep_time = Ts;           
range = 32;                 % Target actual range
velocity = 54;              % m/s
noise = 0;

%% === 2. FMCW Waveform ===
waveform = phased.FMCWWaveform( ...
    'SweepTime', sweep_time, ...
    'SweepBandwidth', bw, ...
    'SampleRate', fs, ...
    'SweepDirection', 'Up', ...           % FIX: no doubled chirps
    'SweepInterval', 'Symmetric', ...
    'NumSweeps', 1024);

tx = waveform();
Nsamples_total = numel(tx);               % FIX: correct number of samples
t = (0:Nsamples_total-1)'/fs;             % FIX: column vector matching tx1d

%% === 3. Target, Radar, and Channel Setup ===
rcs_val = 0.020;
target = phased.RadarTarget('MeanRCS', rcs_val, 'OperatingFrequency', fc);
channel = phased.FreeSpace( ...
    'OperatingFrequency',fc, ...
    'TwoWayPropagation',true, ...
    'SampleRate',fs);

target_motion = phased.Platform('InitialPosition',[range;0;0],'Velocity',[velocity;0;0]);
radar_motion = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);

collector = phased.Collector('OperatingFrequency',fc);
radiator = phased.Radiator('OperatingFrequency',fc);
transmitter = phased.Transmitter('PeakPower',1,'Gain',30);

[tgt_pos, tgt_vel] = target_motion(sweep_time);
[radar_pos, radar_vel] = radar_motion(sweep_time);
[~, ang] = rangeangle(tgt_pos, radar_pos);

%% === 4. Tx → Channel → Rx ===
tx_pwr = transmitter(tx);
tx_radiated = radiator(tx_pwr, ang);
prop_sig = channel(tx_radiated, radar_pos, tgt_pos, radar_vel, tgt_vel);
reflected = target(prop_sig);
rx = collector(reflected, ang);

%% === Add LNA Amplification Stage ===
LNA_gain = 1e6;            % Gain ≈ 120 dB
rx = rx * LNA_gain;

%% === Add noise after amplification ===
rng(0); % Ensures the same noise every run
rx_noisy = awgn(rx, noise, 'measured');

%% === 4.1 Ensure Rx matches Tx in length ===
tx1d = tx(:);                    
rx1d = rx_noisy(:);              

if length(rx1d) > length(tx1d)
    rx1d = rx1d(1:length(tx1d));     % Trim if longer
elseif length(rx1d) < length(tx1d)
    rx1d(end+1:length(tx1d)) = 0;    % Zero-pad if shorter
end

t_rx = (0:length(rx1d)-1)/fs;        % Time vector for Rx

%% === 5. Tx Visualization ===
f_if = 150e6;
tx_if = real(tx1d .* exp(1j*2*pi*f_if*t));

figure;
n_plot = min(2000, length(tx_if));
subplot(5,1,1);
plot(t(1:n_plot)*1e6, tx_if(1:n_plot));
xlabel('Time (\mus)'); ylabel('Amplitude');
title('Tx FMCW Real-Valued IF (Oscilloscope View)');

subplot(5,1,2);
instf_tx = diff(unwrap(angle(tx1d))) * fs / (2*pi);
plot(t(1:end-1)*1e3, instf_tx/1e3);
xlabel('Time (ms)'); ylabel('Freq (kHz)');
title('Tx Instantaneous Frequency');

subplot(5,1,3);
Ytx = abs(fft(tx1d));
faxis = fs*(0:Nsamples_total-1)/Nsamples_total;
plot(faxis/1e3, Ytx); xlim([0 fs/2]/1e3);
xlabel('Freq (kHz)'); ylabel('|FFT|');
title('Tx FFT Spectrum');

win = 256; overlap = 200; nfft = 512;
subplot(5,1,4);
[S_tx, F_tx, T_tx] = spectrogram(tx1d, win, overlap, nfft, fs);
imagesc(T_tx*1e3, F_tx/1e3, mag2db(abs(S_tx)));
axis xy; colorbar;
xlabel('Time (ms)'); ylabel('Freq (kHz)');
title('Tx Spectrogram (Baseband)');

subplot(5,1,5);
[S_rf, F_rf, T_rf] = spectrogram(tx_if, win, overlap, nfft, fs);
imagesc(T_rf*1e6, F_rf/1e6, mag2db(abs(S_rf)));
axis xy; colorbar;
xlabel('Time (\mus)'); ylabel('Freq (MHz)');
title('Tx Spectrogram (Real-Valued IF)');

%% === 6. Rx Visualization (Improved) ===
figure;

% (1) Time-Domain Signal (Normalized & First 2000 Samples)
samples_plot = min(2000, length(rx1d));
subplot(4,1,1);
plot(t_rx(1:samples_plot)*1e6, real(rx1d(1:samples_plot))/max(abs(rx1d)), 'LineWidth',1.2);
xlabel('Time (\mus)'); ylabel('Rx (normalized)');
title('Rx Signal - First 2000 Samples (Time Domain)'); grid on;

% (2) Instantaneous Frequency
instf_rx = gradient(unwrap(angle(rx1d))) * fs / (2*pi);
subplot(4,1,2);
plot(t_rx(1:end-1)*1e3, instf_rx(1:end-1)/1e3, 'LineWidth',1.2);
xlabel('Time (ms)'); ylabel('Freq (kHz)');
title('Rx Instantaneous Frequency'); grid on;

% (3) FFT Spectrum
Yrx = abs(fft(rx1d))/max(abs(rx1d));
faxis_rx = fs*(0:length(rx1d)-1)/length(rx1d);
subplot(4,1,3);
plot(faxis_rx/1e3, Yrx, 'LineWidth',1.2);
xlim([0 fs/2]/1e3); xlabel('Freq (kHz)'); ylabel('|FFT| (norm)');
title('Rx FFT Spectrum'); grid on;

% (4) Spectrogram
subplot(4,1,4);
[S_rx, F_rx, T_rx] = spectrogram(rx1d, win, overlap, nfft, fs);
imagesc(T_rx*1e3, F_rx/1e3, mag2db(abs(S_rx)/max(abs(S_rx(:)))));
axis xy; colorbar;
xlabel('Time (ms)'); ylabel('Freq (kHz)');
title('Rx Spectrogram');

%% === 7. Compute Beat Signal Per Sweep with Doppler Compensation ===
Nsweep = round(sweep_time * fs);
N_sweeps_total = floor(length(tx1d) / Nsweep);

tx_matrix = reshape(tx1d(1:Nsweep * N_sweeps_total), Nsweep, []);
rx_matrix = reshape(rx1d(1:Nsweep * N_sweeps_total), Nsweep, []);

% Perform dechirping per sweep
beat_matrix = zeros(Nsweep, N_sweeps_total);
for k = 1:N_sweeps_total
    beat_matrix(:,k) = rx_matrix(:,k) .* conj(tx_matrix(:,k));
end

% === Estimate Doppler Frequency ===
phase_diff = angle(sum(rx_matrix(:,2:end) .* conj(rx_matrix(:,1:end-1)), 1));
fd_est = mean(phase_diff) / (2*pi*sweep_time);

% === Apply Doppler Phase Compensation ===
doppler_phase = exp(-1j * 2*pi * fd_est * (0:Nsweep-1)' / fs);
beat_corrected = beat_matrix(:,1) .* doppler_phase;

t_beat = (0:Nsweep-1) / fs;

%% === 7.a IF Amplification ===
IF_gain = 1e3;
beat_amp = beat_corrected * IF_gain;

%% === 7.b Bandpass Filter around Expected fr ===
mu = bw / sweep_time;
fr_expected = mu * (2 * range / c);
low_cut = max(1e5, 0.8 * fr_expected);
high_cut = min((fs/2)-1e5, 1.2 * fr_expected);

bpFilt = designfilt('bandpassfir', ...
    'FilterOrder', 80, ...
    'CutoffFrequency1', low_cut, ...
    'CutoffFrequency2', high_cut, ...
    'SampleRate', fs);

beat_filt = filter(bpFilt, beat_amp);

%% === 7.c Plot Beat Signal ===
figure(6);
subplot(3,1,1);
plot(t_beat*1e3, real(beat_filt), 'LineWidth',1.2); grid on;
xlabel('Time (ms)'); ylabel('Re\{Beat\}'); title('Beat Signal - Real');

subplot(3,1,2);
plot(t_beat*1e3, imag(beat_filt), 'LineWidth',1.2); grid on;
xlabel('Time (ms)'); ylabel('Im\{Beat\}'); title('Beat Signal - Imag');

subplot(3,1,3);
plot(t_beat*1e3, abs(beat_filt), 'LineWidth',1.2); grid on;
xlabel('Time (ms)'); ylabel('|Beat|'); title('Beat Signal - Magnitude');

%% === 7.d Instantaneous Frequency - Show Only 10 Chirps ===
num_chirps_to_plot = 15;
samples_to_plot = num_chirps_to_plot * Nsweep;

tx_multi = tx_matrix(:);
rx_multi = rx_matrix(:);

tx_multi = tx_multi(1:samples_to_plot);
rx_multi = rx_multi(1:samples_to_plot);
t_multi = (0:samples_to_plot-1) / fs;

inst_f_tx = gradient(unwrap(angle(tx_multi))) * fs / (2*pi);
inst_f_rx = gradient(unwrap(angle(rx_multi))) * fs / (2*pi);

% تقليل الضوضاء باستخدام moving average بسيط
inst_f_tx = movmean(inst_f_tx, 20);
inst_f_rx = movmean(inst_f_rx, 20);

figure(9);
plot(t_multi*1e6, inst_f_tx/1e6, 'b', 'LineWidth', 1.5); hold on;
plot(t_multi*1e6, inst_f_rx/1e6, 'r', 'LineWidth', 1.2);
xlabel('Time (\mus)'); ylabel('Instantaneous Freq (MHz)');
legend('Tx IF','Rx IF'); grid on;
ylim([min(inst_f_rx/1e6)-20, max(inst_f_tx/1e6)+20]);
title(sprintf('Instantaneous Frequency - First %d Chirps', num_chirps_to_plot));

%% === 8. Range Estimation (Final, Doppler-Compensated) ===
% Windowed FFT on doppler-compensated beat signal
N = length(beat_filt);
beat_win = beat_filt .* hann(N);
Nfft = 4096;

Yb = abs(fft(beat_win, Nfft)) / max(abs(beat_win));
Yb = Yb(1:Nfft/2);
f_axis = fs * (0:(Nfft/2)-1) / Nfft;

% Peak detection with interpolation
[~, idx_peak] = max(Yb);
if idx_peak > 1 && idx_peak < length(Yb)
    alpha = Yb(idx_peak-1); beta = Yb(idx_peak); gamma = Yb(idx_peak+1);
    p = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma);
    fb_est = (idx_peak + p - 1) * (fs / Nfft);
else
    fb_est = f_axis(idx_peak);
end

% === Correct Range Calculation ===
% R = (c * fb * Ts) / (2 * B) → now fb is free from doppler bias
range_est = (c * fb_est * sweep_time) / (2 * bw);

fprintf('Estimated Beat Frequency (No Doppler) = %.2f kHz\n', fb_est/1e3);
fprintf('Estimated Range (Final) = %.2f m\n', range_est);

%% === 9. Doppler Extraction from Phase Progression ===
% Phase difference between consecutive sweeps
phase_diff = angle(sum(rx_matrix(:,2:end) .* conj(rx_matrix(:,1:end-1)), 1));
avg_phase_shift = mean(phase_diff);

% Estimated Doppler frequency and velocity
fd_est = avg_phase_shift / (2*pi*sweep_time);
vel_est_dopp = (fd_est * lambda) / 2;

fprintf('Estimated Doppler Frequency (from phase) = %.2f Hz\n', fd_est);
fprintf('Estimated Velocity (from Doppler) = %.2f m/s\n', vel_est_dopp);


%% === 10. Save Beat Signal to File for Fusion Processing ===
% Choose filename and format
output_filename = 'beat_5_8GHz.mat';  % Change name for each radar
beat_signal = real(beat_filt);        % Save only real part or use abs(beat_filt) if preferred
fs_out = fs;                          % Save sampling rate too
t_out = t_beat;                       % Time vector

% Save as .mat for MATLAB processing
save(output_filename, 'beat_signal', 'fs_out', 't_out');

% (Optional) Also save as .txt for other tools
txt_filename = 'beat_5_8GHz.txt';
writematrix([t_out(:) beat_signal(:)], txt_filename);

fprintf('✔ Beat signal saved to %s and %s\n', output_filename, txt_filename);