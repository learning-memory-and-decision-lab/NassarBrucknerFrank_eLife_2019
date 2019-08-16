%% This script is a tutorial on time frequency analysis copied from here:
% http://bishoptechbits.blogspot.com/2010/07/simplified-introduction-to-time.html

% and extended by MRN to better understand the tradeoffs navigated by
% various wavelets and tapering strategies


%% create simple sin wave.
D = 1; % signal duration, 1 second
S = 1000; % sampling rate, i.e. N points pt sec used to represent sine wave
F = 50; % frequency in Hz
P = .25 % NB. 2 pi radians = 360 degrees
T = 1/S; % sampling period, i.e. for this e.g. points at 1 ms intervals
t = [T:T:D]; % time vector %NB corrected from previous version
myphi=2*pi*P;
mysig = sin(2*F*t*pi+ myphi); % sinusoidal signal; 2*F*t*pi is freq in radians
figure; plot(t,mysig);
xlabel('Time (seconds)');
ylabel('Amplitude');


%% create complex fixed phase sin wave.
D = 1.2; %signal duration
S = 1000; % sampling rate, i.e. N points pt sec used to represent sine wave
F = [200 2 20 40]; % 4 frequencies in Hz
w = 2*pi*F; % convert frequencies to radians
P = [0 .5 .25 .3]; % 4 corresponding phases
A = [.2 .5 .3 .2]; % corresponding amplitudes
T = 1/S; % sampling period, i.e. for this e.g. points at 1 ms intervals
t = [T:T:D]; % time vector %NB this has been corrected from previous version
mysig=zeros(1,length(t)); %initialise mysig
myphi=2*pi*P; % compute phase angle
nfreqs=length(F); % N frequencies in the complex
% Add all sine waves together to give composite
for thisfreq=1:nfreqs
mysig = mysig+A(thisfreq)*(sin(w(thisfreq)*t + myphi(thisfreq)));
end
figure; plot(t,mysig);
xlabel('Time (seconds)');
ylabel('Amplitude');





 %-----------------------------------------------------------------------
% do_FFT.m
%-----------------------------------------------------------------------
% assumes you have already created mysig.
L = length(mysig);
myFFT = fft(mysig,S);
myFFT=myFFT/L; %scale the output to 1
freq = S/2*linspace(0,1,S/2);%create a range with all frequency values
% Plot single-sided amplitude spectrum.
figure; stem(freq,abs(myFFT(1:length(freq))));
xlabel('Frequency (Hz)');
ylabel('Amplitude ');
%-----------------------------------------------------------------------




%% create complex fixed phase sin wave.
D = 2; %signal duration
S = 1000; % sampling rate, i.e. N points pt sec used to represent sine wave
F = [10 20 45]; % 4 frequencies in Hz
w = 2*pi*F; % convert frequencies to radians
P = [0 .0 0]; % 4 corresponding phases
A = [1 1 1]; % corresponding amplitudes
T = 1/S; % sampling period, i.e. for this e.g. points at 1 ms intervals
t = [T:T:D]; % time vector %NB this has been corrected from previous version
mysig=zeros(1,length(t)); %initialise mysig
myphi=2*pi*P; % compute phase angle
nfreqs=length(F); % N frequencies in the complex
% Add all sine waves together to give composite
for thisfreq=1:nfreqs
mysig = mysig+A(thisfreq)*(sin(w(thisfreq)*t + myphi(thisfreq)));
end
figure; plot(t,mysig);
xlabel('Time (seconds)');
ylabel('Amplitude');


%% create a new signal with more 20 hz amplitude and insert it into the old signal
%-----------------------------------------------------------------------
% assumes you have already created mysig
A(2)=2; %amplitude of 2nd frequency increased from 1 to 2
mysig2=zeros(1,length(t)); %initialise mysig2
for thisfreq=1:nfreqs
mysig2 = mysig2+A(thisfreq)*(sin(w(thisfreq)*t + myphi(thisfreq)));
end
mysig(500:1000)=mysig2(500:1000);
t=t-.5; %subtract 500 ms to indicate first 500 ms are baseline, i.e. -ve time
figure; plot(t,mysig);
xlabel('Time (seconds)');
ylabel('Amplitude');

eeglab %load EEGLAB commands if you have not already done so
figure;
[ersp,itc,powbase,times,freqs]=...
newtimef(mysig, 2000,[-500 1500],1000, 0,'plotitc','off');


figure;
imagesc(times, freqs, ersp, [-30 30])
set(gca,'YDir','normal') % y axis is inverted if you don't specify this
h = colorbar;


 figure;[ersp,itc,powbase,times,freqs]=...
newtimef( mysig,2000,[-500 1500],1000, 0,'baseline',NaN,'plotitc','off');


 figure;[ersp,itc,powbase,times,freqs]=...
newtimef( mysig,1000,[-500 500],500, 0,'plotitc','off');


for i=0:4
figure;[ersp,itc,powbase,times,freqs]=...
newtimef( mysig,2000,[-500 1500],1000,i,'plotitc','off','erspmax',30);
end



figure;[ersp,itc,powbase,times,freqs]=...
newtimef( mysig,2000,[-500 1500],1000,[1 5],'plotitc','off','erspmax',30);









 %-----------------------------------------------------------------------
% Make complex signal with bursts of power at 8, 25 and 40 Hz
%-----------------------------------------------------------------------
D = 2; %signal duration 2s
S = 1000; % sampling rate, i.e. N points pt sec used to represent sine wave
F = [8 25 40]; % 4 frequencies in Hz
w = 2*pi*F; % convert frequencies to radians
P = [0 0 0]; % 4 corresponding phases
T = 1/S; % sampling period, i.e. for this e.g. points at 1 ms intervals
t = [T:T:D]; % time vector %corrected from previous version
nfreqs=length(F); % N frequencies in the complex
A=ones(nfreqs,length(t)); %amplitudes initialised at one
A(1,500:700)=3; %boost of increased amplitude at each time pt, nb 500 corresponds to end of baseline period
A(2,900:1100)=3;
A(3,1300:1500)=3;
% Add all sine waves together to give composite
for thispt=1:2000
for thisfreq=1:nfreqs
mysig(thisfreq,thispt) =A(thisfreq,thispt)*sin(w(thisfreq)*t(thispt));
end
end
mysig2=sum(mysig,1);
t=t-.5; %subtract 500 ms to indicate baseline is negative
figure; plot(t,mysig2);
xlabel('Time (seconds)');
ylabel('Amplitude');
for i=0:3
% First plot ERSP with constant values of cycles
scrsz = get(0,'ScreenSize'); %used to adjust figure size; not normally needed
figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4])
[ersp,itc,powbase,times,freqs]=...
newtimef( mysig2,2000,[-500 1500],1000,i,'plotitc','off', 'erspmax',30);
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
i = [1 5];
[ersp,itc,powbase,times,freqs]=...
newtimef( mysig2,2000,[-500 1500],1000,i,'plotitc','off', 'erspmax',30);
%-----------------------------------------------------------------------




 %-----------------------------------------------------------------------------------
% Make complex signal with long and short bursts of power at 8 and 45 Hz
%-----------------------------------------------------------------------------------
D = 2; %signal duration 2s
S = 1000; % sampling rate, i.e. N points pt sec used to represent sine wave
F = [8 8 45 45]; % 4 frequencies in Hz
w = 2*pi*F; % convert frequencies to radians
P = [0 0 0 0]; % 4 corresponding phases
T = 1/S; % sampling period, i.e. for this e.g. points at 1 ms intervals
t = [T:T:D]; % time vector %corrected from previous version
t=t(1:length(t)-1); % takes off extra point arising from starting at zero
nfreqs=length(F); % N frequencies in the complex
A=ones(nfreqs,length(t)); %amplitudes initialised at one
A(1,550:570)=3; %boost of increased amplitude at each time pt
A(2,700:900)=3;
A(3,950:970)=3;
A(4,1100:1300)=3;
% Add all sine waves together to give composite
clear mysig
for thispt=1:length(A)
for thisfreq=1:nfreqs
mysig(thisfreq,thispt) =A(thisfreq,thispt)*sin(w(thisfreq)*t(thispt));
end
end
mysig2=sum(mysig,1);
t=t-.5; %subtract 500 ms to indicate baseline is negative
figure; plot(t,mysig2);
xlabel('Time (seconds)');
ylabel('Amplitude');
scrsz = get(0,'ScreenSize');
for pr=[.5, 1,2,4,8,16] %pad ratio
cyc=[1 ]; % cycles (can experiment with different values)
figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/3])
[ersp,itc,powbase,times,freqs]=...
newtimef( mysig2,length(mysig2),[-500 1500],1000, cyc, 'erspmax', 20, 'plotitc', 'off', 'padratio',pr);
text(52,20,strcat('Pad ratio= ',num2str(pr),': N freqs=',num2str(length(freqs))));
end
%-----------------------------------

