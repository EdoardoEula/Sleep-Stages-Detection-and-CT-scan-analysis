clc
clear
close all
%%
load("record.mat");
EEG = double(record);
clear record
fs_EEG = 512;
N = length(EEG);
EEG = EEG(60*fs_EEG:end);
%% FILTERING

% Removing net interference
[b1,a1] = iirnotch(60/(fs_EEG/2),0.2/(fs_EEG/2)); %this function wants radians!
% figure, freqz(b3,a3);

EEG = filtfilt(b1,a1, EEG);

%Low-pass
filter_order = 2;
f_cut = 30;
wc = f_cut/(fs_EEG/2);
[b1,a1] = butter(filter_order,wc,'low');
%figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b1,a1,EEG);

%High-pass
f_cut = 2;
wc = f_cut/(fs_EEG/2);
[b1,a1] = butter(filter_order,wc,'high');
%figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b1,a1,EEG);

%% Dividing into epochs
% 30 sec
lep = 30;
ep_samples = lep*fs_EEG;
ep = cell(1, floor(length(EEG)/ep_samples));
var = cell(floor(length(EEG)/ep_samples), 4);
perc = cell(floor(length(EEG)/ep_samples), 4);
%var_per = cell(floor(length(EEG)/ep_samples), 4);

%EEG bands
f_vect_high = [14 8 4 0.5];
f_vect_low = [30 13 8 4];
names = ["beta", "alpha", "theta", "delta"];

%% PSD
for i = 1:floor(length(EEG)/ep_samples)
 %for i = 150:152
    ep{i} = EEG(1+ep_samples*(i-1):ep_samples*(i));

    %PSD with fft
        DFTx = fft(ep{i});
        DFTx = DFTx(1:ep_samples/2 + 1);
        PSDx_ep = (1/(fs_EEG*ep_samples)) * abs(DFTx).^2;
        PSDx_ep(2:end) = 2*PSDx_ep(2:end);
        freq1_ep = 0:fs_EEG/ep_samples:fs_EEG/2;

%     %PSD with periodogram
%     nfft = 16384;
%     [pxx_ep,f_ep]=periodogram(ep{i}, rectwin(length(ep{i})), nfft,fs_EEG);


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

%         figure, plot(freq1,PSDx)
%         ylabel("[power/Hz]")
%         title([names{j}, " epoch", i])
%         xlim([0 f_vect_low(j) + 5])
        var{i,j} = trapz(freq1, PSDx);
        perc{i, j} = trapz(freq1, PSDx)/trapz(freq1_ep, PSDx_ep);

%         %PSD with periodogram
%         nfft=16384;
%         [pxx,f]=periodogram(ep_filt, rectwin(length(ep_filt)), nfft,fs_EEG);
%         
% 
% %         figure, 
% %         plot(f,10*log10(pxx))
% %         ylabel("[dB/Hz]")
% %         title([names{j}, " epoch", i])
% %         xlim([0 f_vect_low(j) + 5])
%         var_per{i, j} = trapz(f, pxx)/trapz(f_ep, pxx_ep);
    end
end

var = cell2table(var, 'VariableNames', names);
%save("variances.mat", "var");

perc = cell2table(perc, 'VariableNames', names);
%save("percentages.mat", "perc");

% var_per = cell2table(var_per, 'VariableNames', names);
%save("variances_per.mat", "var_per");

%%
%Plot of the variances for each frequency band
% lim = [0 9];
% figure
% t = (1:lep:lep*floor(length(EEG)/ep_samples))/3600;
% subplot(411), plot(t,var.beta), title("Beta"),ylim([0 1]), xlim(lim);
% subplot(412), plot(t,var.alpha), title("Alpha"),ylim([0 1]), xlim(lim);
% subplot(413), plot(t,var.theta), title("Theta"),ylim([0 1]), xlim(lim);
% subplot(414), plot(t,var.delta), title("Delta"),ylim([0 1]), xlim(lim);
%%
%Plot of the variances for each frequency band
lim = [0 9];
figure
t = (1:lep:lep*floor(length(EEG)/ep_samples))/3600;
subplot(421), plot(t,10*log10(var.beta)), title("Beta"),ylim([-10 50]), xlim(lim);
subplot(423), plot(t,10*log10(var.alpha)), title("Alpha"),ylim([-10 50]), xlim(lim);
subplot(425), plot(t,10*log10(var.theta)), title("Theta"),ylim([-10 50]), xlim(lim);
subplot(427), plot(t,10*log10(var.delta)), title("Delta"),ylim([-10 50]), xlim(lim);

subplot(422), plot(t,perc.beta), title("Beta"),ylim([0 1]), xlim(lim);
subplot(424), plot(t,perc.alpha), title("Alpha"),ylim([0 1]), xlim(lim);
subplot(426), plot(t,perc.theta), title("Theta"),ylim([0 1]), xlim(lim);
subplot(428), plot(t,perc.delta), title("Delta"),ylim([0 1]), xlim(lim);
saveas(gcf(), "Var and perc")