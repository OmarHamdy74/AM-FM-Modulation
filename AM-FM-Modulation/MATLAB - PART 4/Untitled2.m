%%%%%%%%%%%%%%%%%%%%%%%  STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%% 
[y, fs] = audioread('videoplayback (1).wav');

fc = 20000; % carrier frequency
beta1 = 3; 
beta2 = 5;
t = linspace(0, length(y)/fs, length(y));
c = cos(2*pi*fc*t); % Carrier signal

s1 = fmmod(y,fc,fs,beta1);
s2 = fmmod(y,fc,fs,beta2);

% Noise levels to iterate over
noise_levels = [10, 15, 20];

for i = 1:numel(noise_levels)
    snr1 = noise_levels(i);
    snr2 = noise_levels(i);

    s1_noisy = awgn(s1, snr1, 'measured');
    s2_noisy = awgn(s2, snr2, 'measured');

    y1_demod = fmdemod(s1_noisy, fc, fs, beta1);
    y2_demod = fmdemod(s2_noisy, fc, fs, beta2);

    output_snr1 = snr(y, y - y1_demod);
    output_snr2 = snr(y, y - y2_demod);

    fprintf('Output SNR for beta1 = %d, snr = %d: %.2f dB\n', beta1, snr1, output_snr1);
    fprintf('Output SNR for beta2 = %d, snr = %d: %.2f dB\n', beta2, snr2, output_snr2);
end

% Play original and noisy signals
     sound(y, fs);
     pause(length(y)/fs);
     sound(s1_noisy, fs);
     pause(length(s1_noisy)/fs);
     sound(s2_noisy, fs);
     pause(length(s2_noisy)/fs);


     
%%%%%%%%%%%%%%%%%%%%%%%  STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%% 
% Increase deviation ratio
beta1_new = 10; % New deviation ratio for s1
beta2_new = 15; % New deviation ratio for s2

% delta_f1_new = beta1_new * max(abs(y));
% delta_f2_new = beta2_new * max(abs(y));

s1_noisy_high_deviation = fmmod(y, fc, fs, beta1_new);
s2_noisy_high_deviation = fmmod(y, fc, fs, beta2_new);

% FM Demodulation
y1_demod = fmdemod(s1_noisy_high_deviation, fc, fs, beta1_new);
y2_demod = fmdemod(s2_noisy_high_deviation, fc, fs, beta2_new);

% Plot original and demodulated signals for comparison
figure(1);
subplot(3,3,1);
plot(t, y);
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');

figure(1);
subplot(3,3,2);
plot(t, y1_demod);
hold on;
plot(t, y2_demod);
title('Demodulated Signals');
legend('s1 Demodulated', 's2 Demodulated');
xlabel('Time');
ylabel('Amplitude');