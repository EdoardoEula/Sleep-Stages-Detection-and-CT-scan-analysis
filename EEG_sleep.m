clc
clear
close all
%% Loading and cutting the signal
load("record.mat");
EEG = double(record);
clear record
fs_EEG = 512;
EEG = EEG(5*60*fs_EEG:end);
EEG = EEG((length(EEG)/fs_EEG - 5*60)*fs_EEG:-1:1);
N = length(EEG);
%% Filtering
%Removing net interference
[b,a] = iirnotch(60/(fs_EEG/2),0.2/(fs_EEG/2)); %this function wants radians!
% figure, freqz(b3,a3);
EEG = filtfilt(b,a, EEG);

filter_order = 4;

%Low-pass
f_cut = 70;
wc = f_cut/(fs_EEG/2);
[b,a] = butter(filter_order,wc,'low');
%figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b,a,EEG);

%High-pass
f_cut = 0.5;
wc = f_cut/(fs_EEG/2);
[b,a] = butter(filter_order,wc,'high');
%figure; freqz(b1,a1,1024,fs_EEG);

EEG = filtfilt(b,a,EEG);

%plotting power spectrum
noverlap = 20*fs_EEG;
[Pxx, freq1] = pwelch(EEG, hanning(length(EEG)), noverlap, [], fs_EEG);

figure
plot(freq1, Pxx)
xlim([0 15])

%% Dividing into epochs
% 30 sec
lep = 30;
ep_samples = lep*fs_EEG;

%EEG bands
f_vect_high = [8 4 2.5];
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

    %PSD with Welch windowing
    noverlap = 20*fs_EEG;
    [Pxx_ep, freq_ep] = pwelch(ep{i}, hanning(ep_samples), noverlap, [], fs_EEG);
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

        %PSD with Welch windowing
        noverlap = 20*fs_EEG/2;
        [Pxx, freq] = pwelch(ep_filt, hanning(ep_samples), noverlap, [], fs_EEG);

        var{i,j} = trapz(freq, Pxx);
        perc{i, j} = trapz(freq, Pxx)/trapz(freq_ep, Pxx_ep);
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

plot(t,var.alpha),ylim([-10 50]), xlim(lim);
hold on
plot(t,var.theta),ylim([-10 50]), xlim(lim);
plot(t,var.delta),ylim([-10 50]), xlim(lim);

legend('alpha', 'theta', 'delta')

% subplot(322), plot(t,perc.alpha), title("Alpha"),ylim([0 1]), xlim(lim);
% subplot(324), plot(t,perc.theta), title("Theta"),ylim([0 1]), xlim(lim);
% subplot(326), plot(t,perc.delta), title("Delta"),ylim([0 1]), xlim(lim);
% saveas(gcf(), "Var and perc")


%% 
st = 1:5;
domTab = cell(floor(length(EEG)/ep_samples),3);
for i = 1:floor(length(EEG)/ep_samples)
    j = perc{i,:} == max(perc{i,:});
    h = perc{i,:} == min(perc{i,:});
    domTab{i,2} = perc.Properties.VariableNames{j};
    domTab{i,1} = perc.Properties.VariableNames{h};
end



%% hypnogram
stages = [{'awake'}, {'stage 1'}, {'stage 2'}, {'stage 3'}, {'stage 4'}]; 
domArr = zeros(length(domTab), 2);

for i = 1:floor(length(EEG)/ep_samples)
    if domTab{i,1} == "alpha"
        domArr(i,1) = 1;
    end

    if domTab{i,2} == "alpha"
        domArr(i,2) = 1;
    end

    if domTab{i,1} == "theta"
           domArr(i,1) = 2;
    end

    if domTab{i,2} == "theta"
            domArr(i,2) = 2;
    end

    if domTab{i,1} == "delta"
            domArr(i,1) = 3;
    end

    if domTab{i,2} == "delta"
            domArr(i,2) = 3;
    end

end

%alpha 1, theta 2, delta 3

for i = 1:floor(length(EEG)/ep_samples)
    if i+5 < length(domArr)
        a = length(find(domArr(i:i+5, 2) == 1)); %alpha count max
        b = length(find(domArr(i:i+5, 2) == 2)); %theta count max
        c = length(find(domArr(i:i+5, 2) == 3)); %delta count max
        d = length(find(domArr(i:i+5, 1) == 1)); %alpha count min
        e = length(find(domArr(i:i+5, 1) == 3)); %delta count min

        if a >= 3 && e >= 3
            for j = i:i+5
                domTab{i,3} = 4; %awake
            end
        end

        if b >=3 && a >= c
            for j = i:i+5
                domTab{i,3} = 3; %REM
            end
        end

        if c >= 3 && d >= 3
            for j = i:i+5
                domTab{i,3} = 2; %deep sleep
            end
        end

        %         if c > 3
        %             for j = i: i+5
        %                 domTab{j, 3} = 2; %deep sleep
        %             end
        %         end
        %
        %         e = length(find(domArr(i:i+5, 2) == 3));
        %         if e > 2 && b > 3
        %             for j = i: i+5
        %                 domTab{j, 3} = 3; %deep sleep
        %             end
    end

    if isempty(domTab{i, 3})
        domTab{i,3} = 1;
    end

end

%%

for h = 1:1060
    hypno(h) = domTab{h,3};
end

figure
plot(1:1060, hypno)


%%

domTab = cell2table(domTab, 'VariableNames', ["Dominant band", "Stage"]);
save("dom_waves.mat", "domTab");
