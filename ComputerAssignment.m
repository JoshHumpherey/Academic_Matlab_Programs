% Intro to Analog and Digital Communication Systems
% Matlab Computer Assignment
% Written by Josh Humphrey

%% 1: Sampling and Linear Processing
SampleRate = 8; 
SampleInterval = 1/(SampleRate); 
t = -1+SampleInterval:SampleInterval:1-SampleInterval;
G1T = triangularPulse(t);

%1.1: Plot Discrete Version of Signal (Triangle Pulse)
figure(1);
stem(t, G1T)
xlabel('Time (seconds)')
ylabel('G1(t)')
title('Section 1.1')

%1.2: Plot Fourier Spectrum of G1(f) (Triangle Pulse)
Mag1 = abs(fft(G1T));
N = size(t,2);
dF = SampleRate/N;
f = -SampleRate/2:dF:SampleRate/2-dF;

figure(2)
stem(f,Mag1/N);
xlabel('Frequency (Hz)');
ylabel('G1(f)');
title('Signal 1 Magnitude Response');

%1.3: Plot Autocorrelation and Energy Spectral Density (Triangle Pulse)
auto1 = autocorr(G1T);
for i = 1:15
    ESD1(i) = Mag1(i)*Mag1(i);
end

figure(3)
stem(t,auto1);
xlabel('Time (Seconds)');
ylabel('G1(t)');
title('Signal 1 Autocorrelation');

figure(4)
stem(f,ESD1);
xlabel('Frequency (Hz)');
ylabel('G1(f)');
title('Signal 1 ESD');

figure(5)
title('Verification of ESD and Autcorrelation as Fourier Pair');
verify1 = fft(auto1);
subplot(2,1,1);
stem(f,verify1,'r')
xlabel('Frequency (Hz)');
ylabel('G1(f)');
subplot(2,1,2);
stem(f,ESD1,'b')
xlabel('Frequency (Hz)');
ylabel('G1(f)');

%% 1.1: Plot Discrete Version of Signal (Rectangular Pulse)
SampleRate = 8; 
SampleInterval = 1/(SampleRate); 
t = -1+SampleInterval:SampleInterval:1-SampleInterval;
G2T = rectangularPulse(t);

figure(6);
stem(t, G2T)
xlabel('Time (seconds)')
ylabel('G2(t)')
title('Section 1.1')

%1.2: Plot Fourier Spectrum of G1(f) (Rectangular Pulse)
Mag2 = abs(fft(G2T));
N = size(t,2);
dF = SampleRate/N;
f = -SampleRate/2:dF:SampleRate/2-dF;

figure(7)
stem(f,Mag2/N);
xlabel('Frequency (Hz)');
ylabel('G2(f)');
title('Signal 2 Magnitude Response');

%1.3: Plot Autocorrelation and Energy Spectral Density (Rectangular Pulse)
auto2 = autocorr(G2T);
for i = 1:15
    ESD2(i) = Mag2(i)*Mag2(i);
end

figure(8)
stem(t,auto2);
xlabel('Time (Seconds)');
ylabel('G2(t)');
title('Signal 2 Autocorrelation');

figure(9)
stem(f,ESD2);
xlabel('Frequency (Hz)');
ylabel('G2(f)');
title('Signal 2 ESD');

figure(10)
title('Verification of ESD and Autcorrelation as Fourier Pair');
verify1 = fft(auto2);
subplot(2,1,1);
stem(f,verify1,'r')
xlabel('Frequency (Hz)');
ylabel('G2(f)');
subplot(2,1,2);
stem(f,ESD2,'b')
xlabel('Frequency (Hz)');
ylabel('G2(f)');

