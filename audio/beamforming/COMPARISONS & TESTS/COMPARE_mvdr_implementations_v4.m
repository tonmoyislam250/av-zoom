%% ############ FAIR MVDR COMPARISON SCRIPT ############
clear; clc;

addpath('/Users/emonchowdhury/Desktop/Phase 2/av_zoom/audio/beamforming');
addpath('/Users/emonchowdhury/Desktop/Phase 2/av_zoom/audio/beamforming/Taki');
addpath('/Users/emonchowdhury/Desktop/Phase 2/av_zoom/audio/beamforming/Metrics');
addpath('/Users/emonchowdhury/Desktop/Phase 2/av_zoom/audio/beamforming/Test_audio/');

%% ================= CONFIG =================
fs = 16000;
theta_target = 0;
theta_noise  = 40;

theta_target_test = 0;  

SNR_dB       = 0;

mic_pos = [0; 0.08];
c = 343;

clean_path = 'male_clean_15s.wav';
noise_path = 'female_piano_14s.wav';

%% ================= LOAD SIGNALS =================
[s_clean, fs1] = audioread(clean_path);
[v_noise, fs2] = audioread(noise_path);
assert(fs1==fs && fs2==fs);

L = min(length(s_clean), length(v_noise));
s_clean = s_clean(1:L);
v_noise = v_noise(1:L);

%% ================= MULTICHANNEL MIXTURE =================
x = generate_multichannel_test_signal_v2( ...
        clean_path, noise_path, ...
        theta_target, theta_noise, ...
        SNR_dB);

x = x(1:L,:);
x_mono = x(:,1);

%% ================= YOUR MVDR (UNCHANGED) =================
N = 512; hop = 256; nfft = 512;
win = hann(N,'symmetric');
alpha = 0.98;
delta = 1e-3;

micPos = [0; 0.08];
numMics = 2;

sigLen = size(x,1);

x_pad = [zeros(N,numMics); x; zeros(N,numMics)];
X = stft_multichannel(x_pad, win, hop, nfft);

[numFreqs, numFrames, ~] = size(X);
freqs = linspace(0, fs/2, numFreqs);

Rxx = init_covariance(numFreqs, numMics);
for n = 1:numFrames
    Rxx = update_covariance(Rxx, squeeze(X(:,n,:)), alpha);
end

d = compute_steering_vector(theta_target_test, freqs, micPos, c);
W = compute_mvdr_weights(Rxx, d, delta);
Y = apply_mvdr(X, W);

y_mvdr = istft(Y, ...
     'Window', win, ...
     'OverlapLength', length(win)-hop, ...
     'FFTLength', nfft, ...
     'FrequencyRange', 'onesided');

% Crop to original signal length
y_mvdr_yours = real(y_mvdr(1:sigLen));

%% ================= TAKI MVDR (UNCHANGED) =================
Y_taki = mvdr_beamformer(x_pad, fs, theta_target_test, mic_pos);

[~, params] = compute_stft(s_clean, fs);
y_mvdr_taki = compute_istft(Y_taki, params);
% ---------- Crop to original length ----------
y_mvdr_taki = y_mvdr_taki(N+1 : N+sigLen);


%% ================= MATLAB MVDR (REFERENCE) =================
%% ================= MATLAB MVDR (REFERENCE) =================

% ---- MATLAB Subband MVDR Beamformer ----
sa = phased.ULA( ...
        'NumElements', 2, ...
        'ElementSpacing', abs(diff(mic_pos)));

bm = phased.SubbandMVDRBeamformer( ...
    'SensorArray', sa, ...
    'OperatingFrequency', fs/2+1, ...
    'SampleRate', fs, ...
    'Direction', [theta_target_test; 0], ...
    'DiagonalLoadingFactor', delta, ...
    'WeightsOutputPort', false);

% MATLAB expects time-domain multichannel signal
y_mvdr_matlab = bm(x);

% Ensure column vector
y_mvdr_matlab = y_mvdr_matlab(:);

% Trim to match signal length
y_mvdr_matlab = y_mvdr_matlab(1:sigLen);





%% ================= ALIGN SIGNALS =================
Lmin = min([length(s_clean), length(x_mono), ...
            length(y_mvdr_yours), length(y_mvdr_taki), ...
            length(y_mvdr_matlab)]);

s  = s_clean(1:Lmin);
x0 = x_mono(1:Lmin);
y1 = y_mvdr_yours(1:Lmin);
y2 = y_mvdr_taki(1:Lmin);
y3 = real(y_mvdr_matlab(1:Lmin));


%% ================= METRICS =================

% ---- SI-SDR (CORRECT) ----
sisdr_mic   = si_sdr(x0, s);
sisdr_yours = si_sdr(y1, s);
sisdr_taki  = si_sdr(y2, s);
sisdr_matlab = si_sdr(y3, s);


% ---- STOI (FAIR) ----
stoi_mic   = stoi(s, x0, fs);
stoi_yours = stoi(s, y1, fs);
stoi_taki  = stoi(s, y2, fs);
stoi_matlab = stoi(s, y3, fs);

% ================= OSINR =================
% residuals
res_mic   = x0 - s;
res_yours = y1 - s;
res_taki  = y2 - s;
res_matlab = y3 - s;

osinr_mic   = 10*log10( sum(s.^2) / (sum(res_mic.^2)   + 1e-12) );
osinr_yours = 10*log10( sum(s.^2) / (sum(res_yours.^2) + 1e-12) );
osinr_taki  = 10*log10( sum(s.^2) / (sum(res_taki.^2)  + 1e-12) );
osinr_matlab = 10*log10( sum(s.^2) / (sum(res_matlab.^2) + 1e-12) );

%% ================= DISPLAY =================
fprintf('\n===== FAIR MVDR COMPARISON =====\n');

fprintf('\nSI-SDR (dB):\n');
fprintf('Mic:     %.2f\n', sisdr_mic);
fprintf('Yours:   %.2f\n', sisdr_yours);
fprintf('Taki:    %.2f\n', sisdr_taki);
fprintf('MATLAB:  %.2f\n', sisdr_matlab);

fprintf('\nSTOI:\n');
fprintf('Mic:     %.3f\n', stoi_mic);
fprintf('Yours:   %.3f\n', stoi_yours);
fprintf('Taki:    %.3f\n', stoi_taki);
fprintf('MATLAB:  %.3f\n', stoi_matlab);

fprintf('\nOSINR (dB):\n');
fprintf('Mic:     %.2f\n', osinr_mic);
fprintf('Yours:   %.2f\n', osinr_yours);
fprintf('Taki:    %.2f\n', osinr_taki);
fprintf('MATLAB:  %.2f\n', osinr_matlab);






%% $$$$$$$$$$$$$ save outputs $$$$$$$$$$$$$$$$$$$$$$$
audiowrite('../Test_output/COMPARE1_V4_mvdr_mic1_output.wav', s_clean, fs);
audiowrite('../Test_output/COMPARE1_V4_mvdr_emon_output.wav', y_mvdr_yours, fs);
audiowrite('../Test_output/COMPARE1_V4_mvdr_taki_output.wav', y_mvdr_taki, fs);
audiowrite('../Test_output/COMPARE1_V4_mvdr_matlab_output.wav', real(y_mvdr_matlab), fs);
