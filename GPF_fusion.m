%% === Global Power Fusion (GPF – Strict Paper Version) ===
clc; clear; close all;

% === 1. Fixed Fusion Weights (example based on SNR ratios) ===
w1 = 0.1; 
w2 = 0.9;

%% === 2. Load Beat Signals ===
load('beat_5_8GHz.mat','beat_signal','fs_out','t_out');
s1 = beat_signal; fs = fs_out; t = t_out;

load('beat_24GHz.mat','beat_signal');
s2 = beat_signal;

%% === 3. Radar Parameters ===
c = 3e8; B = 150e6; Ts = 667e-9; mu = B/Ts;
fprintf('\n--- GPF (Strict Power-Domain) ---\n');
fprintf('B=%.1f MHz | Ts=%.3e s | mu=%.3e Hz/s | fs=%.1f MHz\n', ...
    B/1e6, Ts, mu, fs/1e6);

%% === 4. Estimate fb for each radar separately ===
fb_5_8 = estimate_fb_simple(s1, fs);
fb_24  = estimate_fb_simple(s2, fs);

%% === 5. Normalize and Convert to Power Domain ===
P1 = abs(s1).^2;    P1_n = P1 / max(P1);
P2 = abs(s2).^2;    P2_n = P2 / max(P2);

%% === 6. Global Power Fusion ===
P_fused = w1 * P1_n + w2 * P2_n;       % Power-domain weighted sum
s_fused = sqrt(P_fused);               % Reconstruct amplitude (envelope)

%% === 7. Fused Beat Frequency and Range ===
fb_fused = w1 * fb_5_8 + w2 * fb_24;   % Weighted beat frequency
R_fused  = (c * fb_fused) / (2 * mu);  % Range estimation

%% === 8. Display Results ===
fprintf('fb_5.8GHz = %.3f MHz\n', fb_5_8/1e6);
fprintf('fb_24GHz  = %.3f MHz\n', fb_24/1e6);
fprintf('fb_fused  = %.3f MHz\n', fb_fused/1e6);
fprintf('R_fused   = %.2f m\n', R_fused);

%% === 9. Plot Signals ===
plot_beats_and_fused(s1, s2, s_fused, t);

%% === Helper Functions ===
function fb = estimate_fb_simple(sig, fs)
    N = length(sig);
    Y = abs(fft(sig .* hann(N), 4096));
    f = fs * (0:(4096/2)-1)/4096;
    [~, idx] = max(Y(1:4096/2));
    fb = f(idx);
end

function plot_beats_and_fused(s1, s2, s_fused, t)
    figure;
    subplot(3,1,1); plot(t, abs(s1)); title('Beat Signal 5.8GHz');
    subplot(3,1,2); plot(t, abs(s2)); title('Beat Signal 24GHz');
    subplot(3,1,3); plot(t, s_fused, 'r'); title('Fused Signal (GPF Power-Domain)');
end
