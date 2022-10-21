clc
clear
close all
%%
load("record.mat");
EEG = double(record);
clear record
fs_EEG = 512;
EEG = EEG(5*60*fs_EEG:end);
EEG = EEG((length(EEG)/fs_EEG - 5*60)*fs_EEG:-1:1);
N = length(EEG);
%% Dividing into epochs
% 30 sec
lep = 30;
ep_samples = lep*fs_EEG;

%EEG bands
f_vect_high = [8 4 2];
f_vect_low = [14 8 4];
names = ["alpha", "theta", "delta"];

%% PSD
ep = cell(1, floor(length(EEG)/ep_samples));
var = cell(floor(length(EEG)/ep_samples), 3);
perc = cell(floor(length(EEG)/ep_samples), 3);

for i = 1:floor(length(EEG)/ep_samples)
    %for i = 150:152
    ep{i} = EEG(1+ep_samples*(i-1):ep_samples*(i));
    ep{i} = detrend(ep{i});

    %Removing net interference
    [b,a] = iirnotch(60/(fs_EEG/2),0.2/(fs_EEG/2)); %this function wants radians!
    % figure, freqz(b3,a3);
    ep{i} = filtfilt(b,a, ep{i});

%     %Low-pass
%     filter_order = 2;
%     f_cut = 15;
%     wc = f_cut/(fs_EEG/2);
%     [b,a] = butter(filter_order,wc,'low');
%     %figure; freqz(b1,a1,1024,fs_EEG);
% 
%     ep{i} = filtfilt(b,a,ep{i});
% 
%     %High-pass
%     f_cut = 2;
%     wc = f_cut/(fs_EEG/2);
%     [b,a] = butter(filter_order,wc,'high');
%     %figure; freqz(b1,a1,1024,fs_EEG);
% 
%     ep{i} = filtfilt(b,a,ep{i});

    %PSD with Welch windowing
    noverlap = 20*fs_EEG/2;
    Pxx_ep = pwelch(ep{i}, hanning(ep_samples), noverlap, [], fs_EEG);
    freq1_ep = 0:fs_EEG/length(Pxx_ep):fs_EEG/2;

    for j = 1 : length(f_vect_low)
        filter_order = 4;
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

        %         %PSD with fft
        %         DFTx = fft(ep_filt);
        %         DFTx = DFTx(1:ep_samples/2 + 1);
        %         PSDx = (1/(fs_EEG*ep_samples)) * abs(DFTx).^2;
        %         PSDx(2:end) = 2*PSDx(2:end);
        %         freq1 = 0:fs_EEG/ep_samples:fs_EEG/2;

        %         figure, plot(freq1,PSDx)
        %         ylabel("[power/Hz]")
        %         title([names{j}, " epoch", i])
        %         xlim([0 f_vect_low(j) + 5])

        %PSD with Welch windowing
        noverlap = 20*fs_EEG/2;
        Pxx = pwelch(ep_filt, hanning(ep_samples), noverlap, [], fs_EEG);
        freq1 = 0:fs_EEG/length(Pxx):fs_EEG/2;

        var{i,j} = trapz(Pxx);
        perc{i, j} = trapz(Pxx)/trapz(Pxx_ep);

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

subplot(321), plot(t,10*log10(var.alpha)), title("Alpha"),ylim([-10 50]), xlim(lim);
subplot(323), plot(t,10*log10(var.theta)), title("Theta"),ylim([-10 50]), xlim(lim);
subplot(325), plot(t,10*log10(var.delta)), title("Delta"),ylim([-10 50]), xlim(lim);

subplot(322), plot(t,perc.alpha), title("Alpha"),ylim([0 1]), xlim(lim);
subplot(324), plot(t,perc.theta), title("Theta"),ylim([0 1]), xlim(lim);
subplot(326), plot(t,perc.delta), title("Delta"),ylim([0 1]), xlim(lim);
saveas(gcf(), "Var and perc")

%% 
domTab = cell(floor(length(EEG)/ep_samples),1);
for i = 1:floor(length(EEG)/ep_samples)
    j = perc{i,:} == max(perc{i,:});
    domTab{i} = perc.Properties.VariableNames{j};
end

domTab = cell2table(domTab, 'VariableNames', "Dominant band");
save("dom_waves.mat", "domTab");
