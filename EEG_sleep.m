clc
clear
close all
%%
load("record.mat");
EEG = double(record);
clear record
fs_EEG = 512;
N = length(EEG);
%% FILTERING

% Removing net interference
[b1,a1] = iirnotch(60/(fs_EEG/2),0.2/(fs_EEG/2)); %this function wants radians!
% figure, freqz(b3,a3);

EEG = filtfilt(b1,a1, EEG);

%Low-pass
filter_order = 4;
f_cut = 30;
wc = f_cut/(fs_EEG/2);
[b1,a1] = butter(filter_order,wc,'low');
% figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b1,a1,EEG);

%High-pass
f_cut = 0.5;
wc = f_cut/(fs_EEG/2);
[b1,a1] = butter(filter_order,wc,'high');
% figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b1,a1,EEG);

%% Dividing into epochs
% 30 sec
ep_samples = 30*fs_EEG;
ep = cell(1, round(length(EEG)/ep_samples));
var = cell(round(length(EEG)/ep_samples), 4);

%EEG bands
f_vect_high = [14 8 4 0.5];
f_vect_low = [30 13 8 4];
names = ["beta", "alpha", "theta", "delta"];

for i = 1:length(EEG)/ep_samples
% for i = 150:152
    ep{i} = EEG(1+ep_samples*(i-1):ep_samples*(i));
    for j = 1 : 4
        %Low-pass
        f_cut = f_vect_low(j);
        wc = f_cut/(fs_EEG/2);
        [b1,a1] = butter(filter_order,wc,'low');
        %figure; freqz(b1,a1,1024,fs_EEG);
        ep_lp = filtfilt(b1,a1,ep{i});

        %High-pass
        f_cut = f_vect_high(j);
        wc = f_cut/(fs_EEG/2);
        [b1,a1] = butter(filter_order,wc,'high');
        %figure; freqz(b1,a1,1024,fs_EEG);
        ep_filt = filtfilt(b1,a1,ep_lp);

        %PSD with fft
        DFTx = fft(ep_filt);
        DFTx = DFTx(1:ep_samples/2 + 1);
        PSDx = (1/(fs_EEG*ep_samples)) * abs(DFTx).^2;
        PSDx(2:end) = 2*PSDx(2:end);
        freq1 = 0:fs_EEG/ep_samples:fs_EEG/2;

%         figure, plot(freq1,10*log10(PSDx))
%         ylabel("[dB/Hz]")
%         title([names{j}, " epoch", i])
%         xlim([0 f_vect_low(j) + 5])
        var{i, j} = trapz(freq1, PSDx);
    end
end

var = cell2table(var, 'VariableNames', names);
save("variances.mat", "var");

%%
%Plot of the variances for each frequency band
t = (1:30:30*round(length(EEG)/ep_samples))/3600;
subplot(411), plot(t,10*log10(var.beta)), title("Beta");
subplot(412), plot(t,10*log10(var.alpha)), title("Alpha");
subplot(413), plot(t,10*log10(var.theta)), title("Theta");
subplot(414), plot(t,10*log10(var.delta)), title("Delta");

