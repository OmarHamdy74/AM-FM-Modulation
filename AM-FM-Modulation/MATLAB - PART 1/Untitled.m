 %%%%%%%%%%%%%%%%%%%%%%%  STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%% 
 [y, fs] = audioread('videoplayback (1).wav');
 
 
 %%%%%%%%%%%%%%%%%%%%%%%  STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
 Y = fft(y);
 L = length(y);
 P1 = abs(Y)/L;
 f = linspace(-fs/2, fs/2, L); % Frequency axis
 
 figure(2);
 subplot(3,3,1)
 plot(f,fftshift(P1))
 title('Single Spectrum of Audio')
 xlabel('Freq (Hz)')
 ylabel('|P1(f)|')
 
 %%%%%%%%%%%%%%%%%%%%%%%  STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%
 duration = length(y)/fs;           %length of the audio file in seconds
 
 %time domain
 t = linspace(0,duration,length(y));
 figure(1);
 subplot(3,3,1)
 plot(t,y);
 xlabel('Time (s)');
 ylabel('Amplitude');
 title('Signal in the Time Domain');
 
 %frequency domain
 Y = fft(y);
 L = length(y);P1 = abs(Y)/L;
 f = linspace(-fs/2, fs/2, L); 
 figure(2);
 subplot(3,3,2)
 plot(f,fftshift(P1))
 title('Signal in the Frequency Domain')
 xlabel('Freq (Hz)')
 ylabel('|P1(f)|')

 
 %%%%%%%%%%%%%%%%%%%%%%%  STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%%%
 fc = 3400; 
 [b,a] = butter(5,fc/(fs/2),'low');     % low-pass filter
 filtered_Data = filter(b,a,y);         
 
 
 sound(y, fs);                % Play the original audio
 pause(length(y)/fs);         
 sound(filtered_Data, fs);    % Play the filtered audio 
 
 % Time Domain
 t = linspace(0, length(y)/fs, length(y));
 figure(1);
 subplot(3,3,2)
 plot(t, filtered_Data)
 title('Filtered Signal in Time Domain')
 xlabel('Time (s)')
 ylabel('Amplitude')
 
% Frequency Domain
 Y_filtered = fft(filtered_Data);
 L = length(filtered_Data);
 P1_filtered = abs(Y_filtered)/L;
 f_filtered = linspace(-fs/2, fs/2, L); 
 
 figure(2);
 subplot(3,3,3)
 plot(f_filtered, fftshift(P1_filtered))
 title('Filtered Signal in Frequency Domain')
 xlabel('Freq (Hz)')
 ylabel('Magnitude')
 

 %%%%%%%%%%%%%%%%%%%%%%%  STEP 5 %%%%%%%%%%%%%%%%%%%%%%%%%%
 noise = 0.1*randn(size(y));   % white noise
 
 % Different cutoff frequencies 
 cutoffs = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000];
 snrs = zeros(size(cutoffs));
 for i = 1:length(cutoffs)
     fc = cutoffs(i); % Cutoff frequency
     [b,a] = butter(5,fc/(fs/2),'low');
     filtered_Data = filter(b,a,y);
     noise_Filtered = filter(b,a,noise);
     snr = 20*log10(norm(filtered_Data)/norm(noise_Filtered));
     snrs(i) = snr;
 end
 
 % Plot the SNRs 
 figure;
 plot(cutoffs, snrs, '-o');
 xlabel('Cutoff Freq (Hz)');
 ylabel('SNR (dB)');
 title('Effect of Cutoff Frequency on SNR Ratio');
 
 
 %%%%%%%%%%%%%%%%%%%%%%%  STEP 6 %%%%%%%%%%%%%%%%%%%%%%%%%%
 pause(length(y)/fs);
 
 fs = 44100;                   % Sampling frequency
 number_of_Bits = 16;          % Number of bits per sample
 number_of_Channels = 1;       % Number of channels (mono)
 
 
 recorder_Object = audiorecorder(fs, number_of_Bits, number_of_Channels);
 
 
 disp('Say "f, s, b, d, n, m"')
 recordblocking(recorder_Object, 6); % Record for 6 seconds
 disp('End of recording.');
 
 audio_Data = getaudiodata(recorder_Object);
 
 sounds = {'f', 's', 'b', 'd', 'n', 'm'};    % test the sounds
 
 % Try different cutoff frequencies and test the sounds
 cutoffs = [1000, 5000, 9000];
 Results = zeros(length(sounds), length(cutoffs));
 for i = 1:length(cutoffs)
     fc = cutoffs(i); % Cutoff frequency
     [b,a] = butter(5,fc/(fs/2),'low');
     filtered_Data = filter(b,a,audio_Data);
     for j = 1:length(sounds)
         sound_Data = filtered_Data((j-1)*fs+1:j*fs);
         soundsc(sound_Data, fs);
         Results(j,i) = input(sprintf('Can you hear the sound "%s" with cutoff %d Hz? (0-10) ', sounds{j}, fc));
     end
 end
 
 % Plot the Results
 figure;
 for j = 1:length(sounds)
     subplot(length(sounds),1,j)
     plot(cutoffs, Results(j,:), '-o');
     xlabel('Cutoff Frequency (Hz)');
     ylabel('Intelligibility (0-10)');
     title(sprintf('Effect of Cutoff Frequency on %s', sounds{j}));
 end
 
  %%%%%%%%%%%%%%%%%%%%%%%  STEP 7 %%%%%%%%%%%%%%%%%%%%%%%%%%
[y, fs] = audioread('videoplayback (1).wav');
%new_fs = 2*fs;

%fc = 48000; % carrier frequency
fc = 20000;
beta = 0.8; % modulation index

%filtered_Data_inter = interp(filtered_Data,new_fs);

s = modulate(y,fc,fs,'amdsb-tc',0.8);
t = linspace(0, length(y)/fs, length(y));
% in time domain
figure(1);
subplot(3,3,3);
plot(t,s);
xlabel('Time (s)')
ylabel('Amplitude')
title('Modulated Signal in Time Domain')

% in frequency domain
S = fft(s);
L = length(s);
P1 = abs(S)/L;
f = linspace(-fs/2, fs/2, L); % Frequency axis
figure(2);
subplot(3,3,4)
plot(f,fftshift(P1));
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Modulated Signal in Frequency Domain')


%%%%%%%%%%%%%%%%%%%%%%%  STEP 8 %%%%%%%%%%%%%%%%%%%%%%%%%%
rectifier_detector = abs(s);             % to rectify the modulated signal
f_cutoff = 500;              %  low-pass filter 
[b, a] = butter(6, f_cutoff/(fs/2), 'low');
filtered_s = filter(b, a, rectifier_detector);

sound(filtered_s,fs);


%%%%%%%%%%%%%%%%%%%%%%%  STEP 9 %%%%%%%%%%%%%%%%%%%%%%%%%%
filtered_s = filtered_s - mean(filtered_s);             % to remove dc component from demodulated signal

E_m = sum(y.^2);
E_s = sum(filtered_s.^2);                      % to compute the energy for 2 signal

k = sqrt(E_m/E_s);                             % to compute the scaling factor

filtered_s_scaled = k * filtered_s;
pause(length(filtered_s)/fs);  
sound(filtered_s_scaled,fs);


%%%%%%%%%%%%%%%%%%%%%%%  STEP 10 %%%%%%%%%%%%%%%%%%%%%%%%%%
% In the time domain
t = linspace(0, length(y)/fs, length(y));
pause(length(filtered_s_scaled)/fs); 
figure(1);
subplot(3,3,4)
plot(t, filtered_s_scaled)
xlabel('Time (s)')
ylabel('Amplitude')
title('Demodulated Signal in Time Domain')

% In the frequency domain
X = fft(filtered_s_scaled);
L = length(filtered_s_scaled);
P1 = abs(X)/L;
f = linspace(-fs/2, fs/2, L);
figure(2);
subplot(3,3,5)
plot(f,fftshift(P1));
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Demodulated Signal in Frequency Domain')