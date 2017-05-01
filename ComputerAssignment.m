%%%%%%%%%%%%%%%%%
%%% APPENDIX %%%%
%%%%%%%%%%%%%%%%%

% Intro to Analog and Digital Communication Systems
% Matlab Computer Assignment
% Written by Josh Humphrey

%% Section I: Sampling and Linear Processing
SampleRate1 = 16; 
SampleInterval1 = 1/SampleRate1;
t = -1+SampleInterval1:SampleInterval1:1-SampleInterval1;
G1T = triangularPulse(t);
legnth = length(t);

%1.1: Plot Discrete Version of Signal (Triangle Pulse)
fig1 = figure(1);
plot(t, G1T)
xlabel('Time (seconds)')
ylabel('G1(t)')
title('Time Domain Traingular Pulse Sample')
print(fig1)

%1.2: Plot Fourier Spectrum of G1(f) (Triangle Pulse) 
trans1 = fft(G1T)*SampleInterval1;
G1F = fftshift(trans1);
f = -1+SampleInterval1:SampleInterval1:1-SampleInterval1;              

fig2 = figure(2)
plot(f,abs(G1F));
xlabel('Frequency (Hz)');
ylabel('G1(f)');
title('Signal 1 Frequency Response');
axis([-0.6 0.6 -0.1 1.1]);
print(fig2)
%1.3: Plot Autocorrelation and Energy Spectral Density (Triangle Pulse)

ESD1=(trans1).*conj(trans1); 
ESD1=fftshift(ESD1);

fig3 = figure(3)
plot(f,ESD1)
xlabel('Frequency (Hz)');
ylabel('Energy');
axis([-0.25 0.25 -0.1 1.1]);
title('ESD of Signal 1');
print(fig3)

temp1 = autocorr(G1T);
for i = 1:21
    auto1(i) = temp1(i);
end

t_scaled = (5:25);
f_scaled = (5:25);

fig4 = figure(4)
plot(t_scaled ,auto1);
xlabel('Frequency (Hz)');
ylabel('G1(f)');
title('Signal 1 Autocorrelation');
print(fig4)

placeholder1 = fft(auto1)*SampleInterval1;
verify1 = fftshift(placeholder1);

fig5 = figure(5)
title('Verification of ESD and Autcorrelation as Fourier Pair');
subplot(2,1,1);
plot(f_scaled,verify1,'r')
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Fourier Transform of Autocorrelation');
subplot(2,1,2);
plot(f,ESD1,'b')
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Energy Spectral Density of Tri Pulse');
print(fig5)

mag = abs(G1T);
Energy = sum(mag.^2)

%% 1.4: Plot Discrete Version of Signal (Rectangular Pulse)
SampleRate1 = 16; 
SampleInterval1 = 1/SampleRate1;
t = -1+SampleInterval1:SampleInterval1:1-SampleInterval1;
legnth = length(t);
G2T = rectangularPulse(t);

fig6 = figure(6);
plot(t, G2T)
xlabel('Time (seconds)'); 
ylabel('G2(t)');
title('Standard Rectangular Pulse');
print(fig6)

%1.5: Plot Fourier Spectrum of G1(f) (Rectangular Pulse)
trans2 = fft(G2T)*SampleInterval1;
G2F = fftshift(trans2);
f = -1+SampleInterval1:SampleInterval1:1-SampleInterval1;

fig7 = figure(7)
plot(f,abs(G2F));
xlabel('Frequency (Hz)');
ylabel('G2(f)');
title('Signal 2 Fourier Transform');
print(fig7)

%1.6: Plot Autocorrelation and Energy Spectral Density (Rectangular Pulse)
ESD2=(trans2).*conj(trans2); 
ESD2=fftshift(ESD2);

fig8 = figure(8)
plot(f,ESD2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Signal 2 ESD');
print(fig8)

temp2 = autocorr(G2T);
for i = 1:21
    auto2(i) = temp2(i);
end

fig9 = figure(9)
plot(t_scaled,auto2);
xlabel('Frequency (Hz)');
ylabel('G2(f)');
title('Signal 2 Autocorrelation');
print(fig9)

placeholder2 = fft(auto2)*SampleInterval1;
verify2 = fftshift(placeholder2);

fig10 = figure(10)
title('Verification of ESD and Autcorrelation as Fourier Pair');
subplot(2,1,1);
plot(f_scaled,verify2,'r')
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Fourier Transform of Autocorrelation');
subplot(2,1,2);
plot(f,ESD2,'b')
xlabel('Frequency (Hz)');
ylabel('G2(f)');
title('Energy Spectral Density of Rect Pulse')
print(fig10)

mag = abs(G2T);
Energy = sum(mag.^2)

%% 1.7: Calculate and Plot Y1(t) and Y2(t)
clc
t=-2:1/16:2;
f = (-8:0.25:8);
Ts = 1/16;
G_1T = triangularPulse(t);            
G_2T = rectangularPulse(t);

int1 = fft(G_1T)*Ts;
int2 = fft(G_2T)*Ts;

G_1F = fftshift(int1);
G_2F = fftshift(int2);

% Low Pass Filter Design
lpf = rectpuls(f,2);
check = 0;
for index = 1:length(f)

    if (lpf(index) == 1)
        check = 1;
    end
    if (lpf(index) == 0) && (check == 1) 
        lpf(index) = 1;
        check = 0;
    end

end

fig11 = figure(11);
plot(f,lpf)
axis([-2 2 -2 2]);
xlabel('Frequency');
ylabel('Amplitude');
title('Filter Response');
print(fig11)

% Calculating and Plotting Y1(f) and Y2(f)
Y_1F = lpf.*G_1F;
Y_2F = lpf.*G_2F;

fig12 = figure(12)
plot(f,Y_1F);
axis([-2 2 -2 2]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Y1(f)');
print(fig12)

fig13 = figure(13)
plot(f,Y_2F);
axis([-2 2 -2 2]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Y2(f)');
print(fig13)

% Calculating and plotting Y1(t) and Y2(t)
Y_1T=ifft(ifftshift(Y_1F));
Y_2T=ifft(ifftshift(Y_2F));

fig14 = figure(14)
plot(t,Y_1T);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Y1(t)');
print(fig14)

fig15 = figure(15)
plot(t,Y_2T);
xlabel('Time (Seconds)');
ylabel('Time (Seconds)');
title('Y2(t)');
print(fig15)

% Calculating and Plotting Y(t) and Y(f)
Y_T = Y_1T.*Y_2T;
fig16 = figure(16);
plot(t,Y_T)
xlabel('Time (Seconds)');
ylabel('Amplitude');
title('Y(t)');
print(fig16)

Y_F = Y_1F.*Y_2F;
fig17 = figure(17);
plot(f,Y_F);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
axis([-1 1 -0.1 1.1]);
title('Y(f)');
print(fig17)

%% Section II: Amplitude Modulation

SampleRate3 = 300
SampleInterval3 = 1/(SampleRate3); 
t = -1+SampleInterval3:SampleInterval3:1-SampleInterval3;
G3T = 2*triangularPulse((t+0.02)/0.02)-3*triangularPulse((t-0.02)/0.02);
n = length(G3T);
dF = SampleRate3/n;
f = -SampleRate3/2:dF:SampleRate3/2-dF;

% Section 2.1: Plotting the Discrete Signal
fig18 = figure(18);
plot(t, G3T);
xlabel('Time (seconds)')
ylabel('G3(t)')
title('Section 2.1 Original Time Domain Signal')
axis([-0.1 0.1 -4 4]);
print(fig18)

n = length(G3T);
nfft = 2^nextpow2(n);
Z = fft(G3T,nfft);
ZZ = Z(1:nfft/2);
nfft3 = SampleRate3*(0:nfft/2-1)/nfft;

fig19 = figure(19)
plot(nfft3,abs(ZZ));
xlabel('Frequency (Hz)');
ylabel('G3(f)');
title('Section 2.1 Original Frequency Domain Signal');
axis([0 25 0 30]);
print(fig19)

mag3 = abs(G3T);
Energy1 = sum(mag3.^2)

mag4 = abs(ZZ);
Energy2 = sum(mag4.^2)

%% Section 2.2: DSB-SC Modulation Section
Carrier_Signal_Amplitude = 1;
Message_Signal_Amplitude = 1;
Fm = 1;
Fc = 300;
m = Message_Signal_Amplitude/Carrier_Signal_Amplitude;

% Representation of the Message Signal
fig20 = figure(20)
Message_Signal = G3T;
subplot (3,1,1);
plot(t,Message_Signal,'b');
axis([-0.1 0.1 -4 4]);
xlabel ('Time (Seconds)');
ylabel ('Amplitude');
title ('Message Signal');

% Representation of the Carrier Signal
Carrier_Signal = sin(2*pi*Fc*t);
subplot (3,1,2);
plot(t,Carrier_Signal,'r');
axis([-0.1 0.1 -(5*10^-13) (5*10^-13)]);
xlabel ('Time (Seconds)');
ylabel ('Amplitude');
title ('Carrier Signal');

% Representation of the DSB-SC Signal
DSBSC_Signal = G3T.*sin(2*pi*Fc*t);
subplot (3,1,3);
plot(t,DSBSC_Signal,'black');
axis([-0.1 0.1 -(5*10^-13) (5*10^-13)]);
xlabel ('Time (Seconds)');
ylabel ('Amplitude');
title ('DSB-SC Signal');
print(fig20)

%% Designing both filters
n = 5;
Wn = 0.5;
[z,p,k] = butter(n,Wn);
trans5 = zp2sos(z,p,k);
lpf5 = freqz(trans5,599);

n = 40;
Wn = 0.5;
[z,p,k] = butter(n,Wn);
trans40 = zp2sos(z,p,k);
lpf40 = freqz(trans40,599);

% Passing the DSB-SC signal through each filter
for i = 1:599
    Result5(i) = lpf5(i).*DSBSC_Signal(i);
end

for i = 1:599
    Result40(i) = lpf40(i).*DSBSC_Signal(i);
end

fig21 = figure(21)
plot(f,Result5);
xlabel('Frequency (Hz)');
ylabel('Magnitude')
title('5th Order Low Pass Filter Result');
print(fig21)

clc

fig22 = figure(22)
plot(f,Result40);
xlabel('Frequency (Hz)');
ylabel('Magnitude')
title('40th Order Low Pass Filter Result');
print(fig22)

% Taking the inverse fourier transform of the response
Result5_TimeDomain =ifft(ifftshift(Result5));
Result40_TimeDomain =ifft(ifftshift(Result40));

fig23 = figure(23)
plot(f,Result5_TimeDomain);
xlabel('Time (Seconds)');
ylabel('Magnitude')
title('5th Order Low Pass Filter Result');
print(fig23)


clc

fig24 = figure(24)
plot(f,Result40);
xlabel('Time (Seconds)');
ylabel('Magnitude')
title('40th Order Low Pass Filter Result');
print(fig24)

%% Section 2.3: DSB+C Modulation
Sig_Amp = 1;
Carrier_Amp = 1;
Fm = 100;
Fc = 300;
m = Sig_Amp/Carrier_Amp;
A = 2;

% Representation of the Message Signal
fig25 = figure(25)
Message_Signal2 = G3T;
subplot (3,1,1);
plot (t,Message_Signal2,'b');
axis([-0.1 0.1 -4 4]);
xlabel ('Time (Seconds)');
ylabel ('Amplitude');
title ('Message Signal');

% Representation of the Carrier Signal
Carrier_Signal = sin(2*pi*Fc*t);
subplot (3,1,2);
plot (t,Carrier_Signal,'r');
axis([-0.1 0.1 -(5*10^-13) (5*10^-13)]);
xlabel ('Time (Seconds)');
ylabel ('Amplitude');
title ('Carrier Signal');

% Representation of the DSB+C Signal
for i = 1:599
    DSBC_Signal(i) = A.*Carrier_Signal(i) + (Carrier_Signal(i).*Message_Signal2(i));
end

subplot (3,1,3);
plot (t,DSBC_Signal,'black');
axis([-0.1 0.1 -(5*10^-13) (5*10^-13)]);
xlabel ('Time (Seconds)');
ylabel ('Amplitude');
title ('DSB+C Signal');
print(fig25)

% Finding the filtered response
n = 40;
Wn = 0.5;
[z,p,k] = butter(n,Wn);
NEW_trans40 = zp2sos(z,p,k);
NEW_lpf40 = freqz(NEW_trans40,599);

for i = 1:599
    NEW_Result40(i) = NEW_lpf40(i).*DSBC_Signal(i);
end

fig26 = figure(26)
plot(f,NEW_Result40);
xlabel('Frequency (Hz)');
ylabel('Magnitude')
title('40th Order Low Pass Filter Result');
axis([-175 175 -4 4*10^-13])
print(fig26)

NEW_Result40_TimeDomain = ifft(ifftshift(NEW_Result40));

fig27 = figure(27)
plot(t,NEW_Result40_TimeDomain);
xlabel('Frequency (Hz)');
ylabel('Magnitude')
title('40th Order Low Pass Filter Result');
print(fig27)