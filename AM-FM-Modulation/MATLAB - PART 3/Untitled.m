%%%%%%%%%%%%%%%%%%%%%%%  STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%% 
[y, fs] = audioread('videoplayback (1).wav');

fc = 20000; % carrier frequency
beta1 = 3; 
beta2 = 5;
t = linspace(0, length(y)/fs, length(y));
c = cos(2*pi*fc*t); % Carrier signal
y = y(1:length(c));

delta_f1 = beta1*max(abs(y)); % Maximum frequency deviation
delta_f2 = beta2*max(abs(y)); % Maximum frequency deviation


s1 = fmmod(y,fc,fs,delta_f1);
s2 = fmmod(y,fc,fs,delta_f2);


% In time domain
figure(1);
subplot(3,3,1);
plot(t,s1);
xlabel('Time (s)')
ylabel('Amplitude')
title('Modulated Signal in Time Domain for Beta = 3')

figure(1);
subplot(3,3,2);
plot(t,s2);
xlabel('Time (s)')
ylabel('Amplitude')
title('Modulated Signal in Time Domain for Beta = 5')


% In frequency domain
S1 = fft(s1);
L = length(s1);
P1 = abs(S1)/L;
f = linspace(-fs/2, fs/2, L); % Frequency axis
figure(2);
subplot(3,3,1)
plot(f,fftshift(P1));
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Modulated Signal in Frequency Domain for Beta = 3')

S2 = fft(s2);
L = length(s2);
P1 = abs(S2)/L;
f = linspace(-fs/2, fs/2, L); % Frequency axis
figure(2);
subplot(3,3,3)
plot(f,fftshift(P1));
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Modulated Signal in Frequency Domain for Beta = 5')



%%%%%%%%%%%%%%%%%%%%%%%  STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%  
local_oscillator = cos(2*pi*fc*t);          % local oscillator signal

% Mix modulated signal and local oscillator signal
mixer_1 = s1 .* local_oscillator;
mixer_2 = s2 .* local_oscillator;

[b,a] = butter(4,2*delta_f1/fs,'low');         % low-pass filter
filtered_1 = filter(b,a,mixer_1);
[b,a] = butter(4,2*delta_f2/fs,'low');         % low-pass filter
filtered_2 = filter(b,a,mixer_2);

% Amplify filtered output to obtain recovered message signal
amplified_1 = 10*filtered_1;
amplified_2 = 10*filtered_2;

% Normalize recovered message signal
y_recov_1 = amplified_1 / max(abs(amplified_1));
sound(y_recov_1, fs);
pause(length(y_recov_1)/fs);
y_recov_2 = amplified_2 / max(abs(amplified_2));
sound(y_recov_2, fs);


%%%%%%%%%%%%%%%%%%%%%%%  STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%  
fs = 44100; 
fc = 20000; 
beta1 = 3;
beta2 = 5; 
t = linspace(0, length(y)/fs, length(y));

%  a sinusoidal tone at 3 kHz
f_tone = 3000; 
y_tone = sin(2*pi*f_tone*t);

% FM signal
delta_f1 = beta1*f_tone; 
y_fm_1 = fmmod(y_tone,fc,fs,delta_f1);

delta_f2 = beta2*f_tone; 
y_fm_2 = fmmod(y_tone,fc,fs,delta_f2);

% plotting of sinusoidal wave
figure(1);
subplot(3,3,3);
plot(t,y_tone);
xlabel('Time (s)');
ylabel('Amplitude');
title('[Sin Wave]Modulating Signal for Beta = 3');

figure(1);
subplot(3,3,4);
plot(t,y_tone);
xlabel('Time (s)');
ylabel('Amplitude');
title('[Sin Wave]Modulating Signal for Beta = 5');


%%%%%%%%%%%%%%%%%%%%%%%  STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%%% 
% in time domain
figure(1);
subplot(3,3,4);
plot(t,y_fm_1);
xlabel('Time (s)');
ylabel('Amplitude');
title('FM Signal for Beta = 3 in time domain');

figure(1);
subplot(3,3,5);
plot(t,y_fm_2);
xlabel('Time (s)');
ylabel('Amplitude');
title('FM Signal for Beta = 5 in time domain');



% In frequency domain
S1 = fft(y_fm_1);
L = length(y_fm_1);
P1 = abs(S1)/L;
f = linspace(-fs/2, fs/2, L); % Frequency axis
figure(2);
subplot(3,3,4)
plot(f,fftshift(P1));
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Modulated Signal in Frequency Domain for Beta = 3 in frequency domain')

S2 = fft(y_fm_2);
L = length(y_fm_2);
P1 = abs(S2)/L;
f = linspace(-fs/2, fs/2, L); % Frequency axis
figure(2);
subplot(3,3,6)
plot(f,fftshift(P1));
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Modulated Signal in Frequency Domain for Beta = 5 in frequency domain')
