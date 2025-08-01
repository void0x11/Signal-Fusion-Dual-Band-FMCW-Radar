%% === Adaptive Power Fusion (APF – Fully Normalized, Eq.24) ===
clc; clear; close all;

%% === 1. Load Beat Signals ===
load('beat_5_8GHz.mat','beat_signal','fs_out','t_out');
s1 = beat_signal; fs = fs_out; t = t_out;

load('beat_24GHz.mat','beat_signal');
s2 = beat_signal;

%% === 2. Parameters ===
c = 3e8; B = 150e6; Ts = 667e-9; mu = B/Ts;

%% === 3. Normalize individual radar powers ===
P1 = abs(s1).^2; P1_n = P1 / max(P1);    % normalized power of radar 1
P2 = abs(s2).^2; P2_n = P2 / max(P2);    % normalized power of radar 2

%% === 4. Parameters for APF ===
segment_len = 256;
N = length(s1);
num_segments = ceil(N / segment_len);

w1_array = zeros(1,num_segments);
w2_array = zeros(1,num_segments);
Wk = zeros(1,num_segments);
fb_seg = zeros(1,num_segments);
P_fused = zeros(N,1);

%% === 5. Estimate individual fb ===
fb_5_8 = estimate_fb_simple(s1, fs);
fb_24  = estimate_fb_simple(s2, fs);

%% === 6. Segment-wise fusion (Eq.22 & Eq.24) ===
for k = 1:num_segments
    idx = (k-1)*segment_len + 1 : min(k*segment_len, N);
    Q1 = sum(P1_n(idx));   % local energy of radar 1
    Q2 = sum(P2_n(idx));   % local energy of radar 2
    
    if (Q1+Q2) > 0
        w1 = Q1 / (Q1 + Q2);
        w2 = Q2 / (Q1 + Q2);
    else
        w1 = 0.5; w2 = 0.5;
    end
    
    w1_array(k) = w1;
    w2_array(k) = w2;
    Wk(k) = Q1 + Q2;                             % segment total energy
    fb_seg(k) = w1 * fb_5_8 + w2 * fb_24;        % segment fused fb
    P_fused(idx) = w1 * P1_n(idx) + w2 * P2_n(idx); % fused power
end

%% === 7. Reconstruct fused signal (normalized) ===
K = max([max(sqrt(P1_n)), max(sqrt(P2_n))]);    % normalization constant
E_fused = K * sqrt(P_fused);                     % Eq.20 & Eq.22

%% === 8. Final fb according to Eq.24 ===
fb_fused = sum(Wk .* fb_seg) / sum(Wk);          % Eq.24
R_fused  = (c * fb_fused) / (2 * mu);

%% === 9. Print results ===
fprintf('\n=== APF (Fully Normalized, Eq.24) ===\n');
fprintf('fb_5.8GHz      = %.3f MHz\n', fb_5_8/1e6);
fprintf('fb_24GHz       = %.3f MHz\n', fb_24/1e6);
fprintf('fb_fused (Eq.24) = %.3f MHz\n', fb_fused/1e6);
fprintf('R_fused        = %.2f m\n', R_fused);

%% === 10. Plot signals ===
figure;
subplot(3,1,1); plot(t, abs(s1)); title('Beat Signal – 5.8 GHz'); grid on;
subplot(3,1,2); plot(t, abs(s2)); title('Beat Signal – 24 GHz'); grid on;
subplot(3,1,3); plot(t, E_fused, 'r'); title('Fused Signal (APF – Normalized)'); grid on;

%% === 11. Plot weights ===
figure;
plot(w1_array,'b-o','LineWidth',1.2); hold on;
plot(w2_array,'r-o','LineWidth',1.2);
xlabel('Segment Index'); ylabel('Weight Value');
legend('w1 (5.8 GHz)','w2 (24 GHz)');
title('APF Segment-wise Weights'); grid on;

%% === Helper Function ===
function fb = estimate_fb_simple(sig, fs)
    N = length(sig);
    Y = abs(fft(sig .* hann(N), 4096));
    f = fs * (0:(4096/2)-1)/4096;
    [~, idx] = max(Y(1:4096/2));
    fb = f(idx);
end