%% === Adaptive Envelope Fusion (AEF – Dynamic fb, بدون Range ثابت) ===
clc; clear; close all;

%% === 1. Load Beat Signals ===
load('beat_5_8GHz.mat','beat_signal','fs_out','t_out');
s1 = beat_signal; fs = fs_out; t = t_out;

load('beat_24GHz.mat','beat_signal');
s2 = beat_signal;

%% === 2. Parameters ===
c = 3e8; B = 150e6; Ts = 667e-9; mu = B / Ts;
fprintf('\n--- Matched Parameters (AEF – Dynamic fb) ---\n');
fprintf('B=%.1f MHz | Ts=%.3e s | mu=%.3e Hz/s | fs=%.1f MHz\n',...
        B/1e6, Ts, mu, fs/1e6);

%% === 3. Estimate fb لكل رادار ديناميكيًا ===
fb_5_8 = estimate_fb_simple(s1, fs);
fb_24  = estimate_fb_simple(s2, fs);

%% === 4. Bandpass Filter لكل رادار باستخدام fb المحسوب ===
bpFilt1 = designfilt('bandpassfir','FilterOrder',80,...
    'CutoffFrequency1',0.8*fb_5_8,'CutoffFrequency2',1.2*fb_5_8,'SampleRate',fs);
bpFilt2 = designfilt('bandpassfir','FilterOrder',80,...
    'CutoffFrequency1',0.8*fb_24,'CutoffFrequency2',1.2*fb_24,'SampleRate',fs);

s1_f = filter(bpFilt1, s1);
s2_f = filter(bpFilt2, s2);

%% === 5. Hilbert + Envelope ===
env1 = abs(hilbert(s1_f));
env2 = abs(hilbert(s2_f));

%% === 6. Normalize Envelopes ===
env1_n = env1 / max(env1);
env2_n = env2 / max(env2);

%% === 7. Adaptive Segment-wise Fusion ===
segment_len = 256;
N = length(env1_n);
num_segments = ceil(N / segment_len);
E_fused = zeros(N,1);
w1_array = zeros(1,num_segments);
w2_array = zeros(1,num_segments);

for k = 1:num_segments
    idx = (k-1)*segment_len + 1 : min(k*segment_len, N);
    seg1 = env1_n(idx);
    seg2 = env2_n(idx);

    Q1 = sum(seg1.^2);
    Q2 = sum(seg2.^2);

    if (Q1 + Q2) > 0
        w1 = Q1 / (Q1 + Q2);
        w2 = Q2 / (Q1 + Q2);
    else
        w1 = 0.5; w2 = 0.5;
    end

    w1_array(k) = w1;
    w2_array(k) = w2;
    E_fused(idx) = w1 * seg1 + w2 * seg2;
end

% Scale back to amplitude
E_fused = E_fused * max([max(env1), max(env2)]);
fprintf('Max amplitude of fused signal: %.3f\n', max(E_fused));

%% === 8. Final Bandpass Filter using متوسط fb ===
fb_center = mean([fb_5_8, fb_24]);
bpFiltF = designfilt('bandpassfir','FilterOrder',80,...
    'CutoffFrequency1',0.8*fb_center,'CutoffFrequency2',1.2*fb_center,'SampleRate',fs);
E_fused_filt = filter(bpFiltF, E_fused);

%% === 9. FFT Processing ===
Nfft = 4096;
win = hann(length(E_fused_filt));
Yb = abs(fft(E_fused_filt .* win, Nfft));
Yb = Yb(1:Nfft/2);
f_axis = fs * (0:(Nfft/2)-1) / Nfft;

%% === 10. Peak Detection & Range Estimation ===
[~, idx_peak] = max(Yb);
fb_fused = f_axis(idx_peak);
R_fused = (c * fb_fused) / (2 * mu);
fprintf('fb_5.8GHz = %.3f MHz\n', fb_5_8/1e6);
fprintf('fb_24GHz  = %.3f MHz\n', fb_24/1e6);
fprintf('fb_fused  = %.3f MHz\n', fb_fused/1e6);
fprintf('R_fused   = %.2f m\n', R_fused);

%% === 11. Plot Envelopes ===
figure;
subplot(3,1,1); plot(t, env1); title('Envelope – Radar 5.8 GHz');
subplot(3,1,2); plot(t, env2); title('Envelope – Radar 24 GHz');
subplot(3,1,3); plot(t, E_fused, 'r'); title('Adaptive Fused Envelope (AEF)');

%% === 12. Plot FFT ===
figure;
plot(f_axis, Yb, 'k','LineWidth',1.5); grid on;
xlabel('Frequency (Hz)'); ylabel('|FFT|');
title('FFT of Adaptive Envelope Fusion Signal');

%% === 13. Plot Segment Weights ===
figure;
plot(w1_array,'b-o','LineWidth',1.2); hold on;
plot(w2_array,'r-o','LineWidth',1.2);
xlabel('Segment Index'); ylabel('Weight Value');
legend('w1 (5.8 GHz)','w2 (24 GHz)');
title('Segment-wise Adaptive Weights');
grid on;

%% === Helper Function ===
function fb = estimate_fb_simple(sig, fs)
    N = length(sig);
    Y = abs(fft(sig .* hann(N), 4096));
    f = fs * (0:(4096/2)-1)/4096;
    [~, idx] = max(Y(1:4096/2));
    fb = f(idx);
end