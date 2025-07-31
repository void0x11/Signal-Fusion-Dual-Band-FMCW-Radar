clc; clear; close all;

%% Radar System Parameters
fc = 5.8e9;               % Carrier Frequency (Hz)
c = 3e8;                  % Speed of light (m/s)
lambda = c/fc;            % Wavelength (m)
range_max = 500;          % Maximum target range (m)
range_res = 1;            % Range resolution (m)

% Compute Sweep Time (Factor of 5.5 for full wave return)
tm = 5.5 * range2time(range_max, c);

% Compute Sweep Bandwidth
bw = rangeres2bw(range_res, c);
sweep_slope = bw / tm;    % Sweep slope (Hz/s)

% Compute Maximum Beat Frequency
fr_max = range2beat(range_max, sweep_slope, c);

% Assume a maximum drone speed of 50 m/s
v_max = 50;
fd_max = speed2dop(2*v_max, lambda); 
fb_max = fr_max + fd_max;

% Sampling Rate: Use max(2*max beat frequency, bandwidth)
fs = max(2*fb_max, bw);

%% Define FMCW Waveform
waveform = phased.FMCWWaveform('SweepTime', tm, ...
    'SweepBandwidth', bw, 'SampleRate', fs);

% Visualize Signal
sig = waveform();
subplot(2,1,1);
plot(0:1/fs:tm-1/fs, real(sig));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW Signal'); axis tight;

subplot(2,1,2);
spectrogram(sig, 32, 16, 32, fs, 'yaxis');
title('FMCW Signal Spectrogram');

%% Define Drone Target Model
drone_range = 50;        % Drone is at 300 m
drone_velocity = 200;      % Drone moving at 20 m/s

drone_rcs = db2pow(10);   % Assumed RCS of a drone (10 dBsm)
drone_target = phased.RadarTarget('MeanRCS', drone_rcs, ...
    'PropagationSpeed', c, 'OperatingFrequency', fc);

% Drone Movement
drone_motion = phased.Platform('InitialPosition', [drone_range;0;50], ...
    'Velocity', [drone_velocity; 0; 0]);

%% Define Propagation Channel (Two-Way Free Space)
channel = phased.FreeSpace('PropagationSpeed', c, ...
    'OperatingFrequency', fc, 'SampleRate', fs, 'TwoWayPropagation', true);

%% Radar System Components
ant_aperture = 0.03;                          % Approximate drone radar antenna aperture (sq. meters)
ant_gain = aperture2gain(ant_aperture, lambda); % Compute gain

tx_power = db2pow(1) * 1e-3;                   % 1 dBm peak power
tx_gain = 20 + ant_gain;                       % Transmitter gain (dB)
rx_gain = 20 + ant_gain;                       % Receiver gain (dB)
rx_nf = 4.5;                                   % Noise figure (dB)

transmitter = phased.Transmitter('PeakPower', tx_power, 'Gain', tx_gain);
receiver = phased.ReceiverPreamp('Gain', rx_gain, 'NoiseFigure', rx_nf, 'SampleRate', fs);

%% Simulate Radar Returns for Multiple Sweeps
Nsweep = 64;
xr = complex(zeros(waveform.SampleRate * waveform.SweepTime, Nsweep));

for m = 1:Nsweep
    % Update drone and radar positions
    radar_motion = phased.Platform('InitialPosition', [0;0;10], 'Velocity', [0;0;0]);
    [radar_pos, radar_vel] = radar_motion(waveform.SweepTime);
    
    % Update drone position (ADD THIS LINE)
    [tgt_pos, tgt_vel] = drone_motion(waveform.SweepTime);  % <-- Fix here

    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);

    % Propagate signal to and from target
    radar_pos = radar_pos(:);
    tgt_pos = tgt_pos(:);  % Now tgt_pos is defined
    radar_vel = radar_vel(:);
    tgt_vel = tgt_vel(:);  % Now tgt_vel is defined
    
    txsig = channel(txsig, radar_pos, tgt_pos, radar_vel, tgt_vel);
    txsig = drone_target(txsig);

    % Receive and process return signal
    txsig = receiver(txsig);
    dechirpsig = dechirp(txsig, sig);

    xr(:,m) = dechirpsig;
end

%% Range-Doppler Processing
rngdopresp = phased.RangeDopplerResponse('PropagationSpeed', c, ...
    'DopplerOutput', 'Speed', 'OperatingFrequency', fc, ...
    'SampleRate', fs, 'RangeMethod', 'FFT', 'SweepSlope', sweep_slope, ...
    'RangeFFTLengthSource', 'Property', 'RangeFFTLength', 2048, ...
    'DopplerFFTLengthSource', 'Property', 'DopplerFFTLength', 256);

figure;
plotResponse(rngdopresp, xr);
xlabel('Speed (m/s)'); ylabel('Range (m)');
title('Range-Doppler Response');

%% CFAR Detection
% Compute Range FFT
range_fft = fft(xr, 2048, 1); % 2048-point FFT
range_spectrum = abs(range_fft).^2; % Compute power spectrum

% Average the range spectrum across sweeps
avg_spectrum = mean(range_spectrum, 2); % Average across columns (sweeps)

% Ensure avg_spectrum is a column vector
avg_spectrum = avg_spectrum(:); % Convert to column vector

% Debug: Check avg_spectrum
disp('Size of avg_spectrum:'); disp(size(avg_spectrum)); % Should be [2048, 1]
disp('Min of avg_spectrum:'); disp(min(avg_spectrum)); % Should be >= 0
disp('Max of avg_spectrum:'); disp(max(avg_spectrum)); % Should be finite

% Configure CFAR Detector
cfar_detector = phased.CFARDetector('NumTrainingCells', 10, 'NumGuardCells', 4, ...
    'ThresholdFactor', 'Auto', 'Method', 'CA', 'OutputFormat', 'Detection index');

% Define probability of false alarm
Pfa = 1e-6; % Probability of false alarm

% Debug: Check Pfa
disp('Pfa:'); disp(Pfa); % Should be a scalar between 0 and 1

% Apply CFAR detector
try
    % Apply CFAR detector to the averaged spectrum
    detections = cfar_detector(avg_spectrum, Pfa); % Only 2 inputs: data and Pfa
catch ME
    % Display error message if CFAR detector fails
    disp('Error in CFAR Detector:');
    disp(ME.message);
    return; % Exit the function if an error occurs
end

% Convert detections to ranges
range_axis = (0:2047) * c / (2*bw); % Range axis (m)
detected_ranges = range_axis(detections); % Convert detection indices to ranges

% Plot results
figure;
plot(range_axis, avg_spectrum); hold on; % Plot the range spectrum
plot(detected_ranges, avg_spectrum(detections), 'ro'); % Plot detected ranges
xlabel('Range (m)'); ylabel('Power'); % Label axes
title('CFAR Detection Results'); % Add title
legend('Range Spectrum', 'Detections'); % Add legend

%% Range Doppler Coupling Effect
fd = -rootmusic(xr(val2ind(rng_est, c/(fs*2)),:), 1, 1/tm);
v_est = dop2speed(fd, lambda) / 2;
deltaR = rdcoupling(fd, sweep_slope, c);
disp(['Estimated Range Error due to Doppler: ', num2str(deltaR), ' m']);

%% Two-Ray Multipath Simulation
txchannel = twoRayChannel('PropagationSpeed', c, 'OperatingFrequency', fc, 'SampleRate', fs);
rxchannel = twoRayChannel('PropagationSpeed', c, 'OperatingFrequency', fc, 'SampleRate', fs);

xr_multipath = helperFMCWTwoRaySimulate(Nsweep, waveform, ...
    phased.Platform('InitialPosition', [0;0;10], 'Velocity', [0;0;0]), ...
    drone_motion, transmitter, txchannel, rxchannel, drone_target, receiver);

% Plot Multipath Response
figure;
plotResponse(rngdopresp, xr_multipath);
xlabel('Speed (m/s)'); ylabel('Range (m)');
title('Range-Doppler Map with Multipath');