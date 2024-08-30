%%%%%%%%%%%%%%%%%%%%%%%  STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%% 
[y, fs] = audioread('videoplayback (1).wav');

fc = 20000; % carrier frequency

s = modulate(y,fc,fs,'amdsb-sc');
% In time domain
t = linspace(0, length(y)/fs, length(y));
figure(1);
subplot(3,3,1);
plot(t,s);
xlabel('Time (s)')
ylabel('Amplitude')
title('Modulated Signal in Time Domain')

% In frequency domain
S = fft(s);
L = length(s);
P1 = abs(S)/L;
f = linspace(-fs/2, fs/2, L); % Frequency axis
figure(2);
subplot(3,3,1)
plot(f,fftshift(P1));
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Modulated Signal in Frequency Domain')

 
 

%%%%%%%%%%%%%%%%%%%%%%%  STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%% 
t = linspace(0, length(s)/fs, length(s));
fc = 20000;      % carrier frequency
phi = 0;         % carrier phase
c = cos(2*pi*fc*t + phi);
y_shifted = s .* c';

f_cutoff = 2000; % cutoff frequency
b = fir1(100, f_cutoff/(fs/2));
message_of_y = filter(b, 1, y_shifted);
sound(y, fs);
pause(length(y)/fs);  
sound(message_of_y, fs);



%%%%%%%%%%%%%%%%%%%%%%%  STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%% 
% In the time domain
t = linspace(0, length(y)/fs, length(y));
figure(1);
subplot(3,3,2)
plot(t, message_of_y)
xlabel('Time (s)')
ylabel('Amplitude')
title('Demodulated Signal in Time Domain')

% In the frequency domain
X = fft(message_of_y);
L = length(message_of_y);
P1 = abs(X)/L;
f = linspace(-fs/2, fs/2, L);
figure(2);
subplot(3,3,2)
plot(f,fftshift(P1));
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Demodulated Signal in Frequency Domain')



%%%%%%%%%%%%%%%%%%%%%%%  STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%%% 
  % Set the frequency offset values to simulate
  some_diiferent_freqs = [0, 100, 1000, 2000];  
  for i = 1:length(some_diiferent_freqs)
      
      freq_offset = some_diiferent_freqs(i);
      carrier_signal = cos(2*pi*(fc+freq_offset)*t);
  
      demodulated_signal = s .* carrier_signal';
  
      filtered_demodulated_signal = filter(b, 1, demodulated_signal); % low-pass filter
  
      % Normalize the filtered demodulated signal
      max_amplitude = max(abs(filtered_demodulated_signal));
      filtered_demodulated_signal = filtered_demodulated_signal / max_amplitude;
  
      pause(length(s)/fs);  
      sound(filtered_demodulated_signal, fs);
  end



%%%%%%%%%%%%%%%%%%%%%%%  STEP 5 %%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the low-pass filter
f_cutoff = 3400; 
b = fir1(100, f_cutoff/(fs/2));

bandlimited_signal = filter(b, 1, y);

max_amplitude = max(abs(bandlimited_signal));
bandlimited_signal = bandlimited_signal / max_amplitude;

fc = 20000; 

t = linspace(0, length(bandlimited_signal)/fs, length(bandlimited_signal));
carrier_signal_1 = exp(1j*2*pi*fc*t);              % for LSB
%carrier_signal_2 = exp(1j*2*pi*f_cutoff*t);          % for USB

% Use the Hilbert transform to obtain the analytic signal of the bandlimited signal
My_signal = hilbert(bandlimited_signal);

% SSB-SC modulation
s_LSB = real(My_signal .* carrier_signal_1');
%s_USB = real(My_signal .* carrier_signal_1' .* carrier_signal_2');

% in the time domain
% for LSB
figure(1);
subplot(3,3,3);
plot(t, s_LSB);
xlabel('Time (s)');
ylabel('Amplitude');
title('Modulated Signal in Time Domain (LSB)');
% % for USB
% figure(1);
% subplot(3,3,4);
% plot(t, s_USB);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Modulated Signal in Time Domain (USB)');

% in the frequency domain
% for LSB
figure(2);
L = length(s_LSB);
f = linspace(-fs/2, fs/2, L); % Frequency axis
S_LSB = fft(s_LSB)/L;           
P_LSB = abs(S_LSB);          
subplot(3,3,3);
plot(f, fftshift(P_LSB));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Modulated Signal in Frequency Domain (LSB)');
% % for USB
% figure(2);
% L = length(s_USB);
% f = linspace(-fs/2, fs/2, L); % Frequency axis
% S_USB = fft(s_USB)/L;           
% P_USB = abs(S_USB);          
% subplot(3,3,4);
% plot(f, fftshift(P_USB));
% xlabel('Frequency (kHz)');
% ylabel('Magnitude');
% title('Modulated Signal in Frequency Domain (USB)');



%%%%%%%%%%%%%%%%%%%%%%%  STEP 6 %%%%%%%%%%%%%%%%%%%%%%%%%%
t = linspace(0, length(s_LSB)/fs, length(s_LSB));
fc = 20000;      % carrier frequency
phi = 0;         % carrier phase
c = cos(2*pi*fc*t + phi);
y_shifted = s_LSB .* c';

f_cutoff = 2000; % cutoff frequency
b = fir1(100, f_cutoff/(fs/2));
message_of_y = filter(b, 1, y_shifted);

pause(length(s)/fs);
sound(message_of_y, fs);


%%%%%%%%%%%%%%%%%%%%%%%  STEP 7 %%%%%%%%%%%%%%%%%%%%%%%%%%
% In the time domain
t = linspace(0, length(s_LSB)/fs, length(s_LSB));
figure(1);
subplot(3,3,5)
plot(t, message_of_y)
xlabel('Time (s)')
ylabel('Amplitude')
title('Demodulated Signal in Time Domain')

% In the frequency domain
X = fft(message_of_y);
L = length(message_of_y);
P1 = abs(X)/L;
f = linspace(-fs/2, fs/2, L);
figure(2);
subplot(3,3,5)
plot(f,fftshift(P1));
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Demodulated Signal in Frequency Domain')




%%%%%%%%%%%%%%%%%%%%%%%  STEP 8 %%%%%%%%%%%%%%%%%%%%%%%%%% 
  % Set the frequency offset values to simulate
  some_diiferent_freqs = [0, 100, 1000, 2000];  
  for i = 1:length(some_diiferent_freqs)
      
      freq_offset = some_diiferent_freqs(i);
      carrier_signal = cos(2*pi*(fc+freq_offset)*t);
  
      demodulated_signal = s_LSB .* carrier_signal';
  
      filtered_demodulated_signal = filter(b, 1, demodulated_signal); % low-pass filter
  
      % Normalize the filtered demodulated signal
      max_amplitude = max(abs(filtered_demodulated_signal));
      filtered_demodulated_signal = filtered_demodulated_signal / max_amplitude;
  
      pause(length(s)/fs);  
      sound(filtered_demodulated_signal, fs);
  end