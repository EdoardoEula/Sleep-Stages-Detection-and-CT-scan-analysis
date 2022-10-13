%EEG
clc
clear
close all
%%%
load("record.mat");
EEG = double(record);
clear record
fs_EEG = 512;

t = 0:1:300;
%% PSD with fft
N = length(EEG); %samples
DFTx = fft(EEG);
DFTx = DFTx(1:N/2+1); 
PSDx = (1/(fs_EEG*N)) * abs(DFTx).^2;
PSDx(2:end) = 2*PSDx(2:end);
freq1 = 0:fs_EEG/N:fs_EEG/2;

figure
plot(freq1,PSDx)

trapz(freq1,PSDx) %check that variance is equal to var(x)

%%
[r, lag] = xcorr(EEG,EEG,'biased');
%figure, stem(lag,r)
nfft = 2^24; %zero-padding
DFTx_indirect = abs(fft(r));
PSDx_indirect = (1/fs_EEG)*DFTx_indirect(1:N/2+1);
PSDx_indirect(2:end) = 2*PSDx_indirect(2:end);
N = length(DFTx_indirect);
freq2 = 0:fs_EEG/N:fs_EEG/2;
[pxx2,f2] =  periodogram(EEG,hamming(length(EEG)),nfft,fs_EEG);

figure
plot(f2,pxx2)

trapz(freq2,PSDx_indirect) %check that variance is equal to var(x)
%%
% Removing net interference
[b1,a1] = iirnotch(60/(fs_EEG/2),0.2/(fs_EEG/2)); %this function wants radians!
% figure, freqz(b3,a3);

EEG = filtfilt(b1,a1, EEG);

% figure
% ha(1) = subplot(211);
% plot(t, EEG(1:1:301));
% xlabel('time [s]'),ylabel('[mV]'),title('EEG')
% 
% ha(2)=subplot(212);
% plot(t, EEG_filtered(1:1:301));
% xlabel('time [s]'),ylabel('[mV]'),title('Cleaned EEG')

% linkaxes(ha,'x')

%% Filtering
%Low-pass
filter_order = 4;
f_cut = 30;
wc = f_cut/(fs_EEG/2);
[b1,a1] = butter(filter_order,wc,'low');
% figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b1,a1,EEG);

% figure
% ha(1) = subplot(211);
% plot(ha(1), t, EEG(5400:1:5700))
% title(ha(1), "Noisy EEG")
% ha(2) = subplot(212);
% plot(ha(2), t, EEG_lp(5400:1:5700))
% title(ha(2), "Low-passed EEG")

%%
%High-pass
f_cut = 0.5;
wc = f_cut/(fs_EEG/2);
[b1,a1] = butter(filter_order,wc,'high');
% figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b1,a1,EEG);

% figure
% ha(1) = subplot(211);
% plot(ha(1), t, EEG(1:1:301))
% title(ha(1), "Noisy EEG")
% ha(2) = subplot(212);
% plot(ha(2), t, EEG_filtered(1:1:301))
% title(ha(2), "Filtered EEG")

%% Frequency spectrum filtered signal

EEG_spectrum = abs(fft(EEG));
f = linspace(0, fs_EEG, length(EEG_spectrum));

figure
plot(f, EEG_spectrum);
xlim([0 fs_EEG/2]);
title("Cleaned EEG")

%% PSD
N = length(EEG); %samples
DFTx = fft(EEG);
DFTx = DFTx(1:N/2+1); 
PSDx = (1/(fs_EEG*N)) * abs(DFTx).^2;
PSDx(2:end) = 2*PSDx(2:end);
freq1 = 0:fs_EEG/N:fs_EEG/2;

figure
plot(freq1,PSDx)

trapz(freq1,PSDx) %check that variance is equal to var(x)

%% Epochs
% 10 min
ep_samples = (10*60)/0.002;
n = 9;
ep0 = cell(n,1);
ep_f = cell(n,1);
for i = 1:n
    ep0{i} = EEG(1+ep_samples*(i-1):ep_samples*(i));
    ep_f{i} = EEG(1+length(EEG)-ep_samples*(n-(i-1)):length(EEG)-ep_samples*(n-i));
end

% ep_mid = EEG((length(EEG)-1)/2:(length(EEG)-1)/2 + dur -1);
% ep_f = EEG(length(EEG)-dur+1:length(EEG));

%% Plot
t = 0.002:0.002:((ep_samples)*0.002);
figure
plot(t, ep_f{1});
%% EEG bands

f_vect_high = [14 8 4 0.5];
f_vect_low = [30 13 8 4];
names = [{"alpha"}, {"theta"}, {"delta"}];
filter_order = 4;
% ep = [{ep0}, {ep_mid}, {ep_f}];

i = 2;
f_lim = fs_EEG/2;
ep = ep0;

for j = 1:n
        %Low-pass
        f_cut = f_vect_low(i);
        wc = f_cut/(fs_EEG/2);
        [b1,a1] = butter(filter_order,wc,'low');
        %figure; freqz(b1,a1,1024,fs_EEG);
        ep_lp = filtfilt(b1,a1,ep{j});

        %High-pass
        f_cut = f_vect_high(i);
        wc = f_cut/(fs_EEG/2);
        [b1,a1] = butter(filter_order,wc,'high');
        %figure; freqz(b1,a1,1024,fs_EEG);
        ep_filt = filtfilt(b1,a1,ep_lp);

        %PSD
        N = length(ep_filt); %samples
        DFTx = fft(ep_filt);
        DFTx = DFTx(1:N/2+1);
        PSDx = (1/(fs_EEG*N)) * abs(DFTx).^2;
        PSDx(2:end) = 2*PSDx(2:end);
        freq1 = 0:fs_EEG/N:fs_EEG/2;

        figure
        plot(freq1,PSDx)
        title(j)
% 
%         t = 0.002:0.002:((ep_samples)*0.002);
%         figure
%         plot(t, ep_filt);

end

%% Plot

%% Epilepsy

for j = 1:n
       %pass band
       f_epil = [2.5, 4];
        
       ep_filt=bandpass(ep{j}, f_epil, fs_EEG);
       %PSD
        N = length(ep_filt); %samples
        DFTx = fft(ep_filt);
        DFTx = DFTx(1:N/2+1);
        PSDx = (1/(fs_EEG*N)) * abs(DFTx).^2;
        PSDx(2:end) = 2*PSDx(2:end);
        freq1 = 0:fs_EEG/N:fs_EEG/2;
        figure
        plot(freq1,PSDx)
        title(j)
end