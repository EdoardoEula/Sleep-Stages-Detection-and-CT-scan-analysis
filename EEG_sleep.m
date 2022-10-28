clc
clear
close all
%% Loading and cutting the signal
load("record.mat");
EEG = double(record);
clear record
fs_EEG = 512;
EEG = EEG(2*60*fs_EEG:end);
EEG = EEG((length(EEG)/fs_EEG - 1*60)*fs_EEG:-1:1);
N = length(EEG);
%% Complete Power Spectrum
noverlap = 0.5*N;
[Pxx, freq1] = pwelch(EEG, hanning(length(EEG)), noverlap, [], fs_EEG);

figure
plot(freq1, Pxx)
xlabel('[Hz]', fontsize=14);
ylabel('[Power/f]',fontsize=14);

%% Filtering
%Removing net interference
[b,a] = iirnotch(60/(fs_EEG/2),0.2/(fs_EEG/2)); %this function wants radians!
figure 
%freqz(b,a);
[hz,hp,hl]=zplane(b,a);

set(hz, markerSize=25, color='b',linewidth=2);
set(hp, markerSize=25, color='r',linewidth=2);
set(hl,linewidth=2)
title('Notch Filter (60 Hz)','FontSize',16)

EEG = filtfilt(b,a, EEG);
%%
filter_order = 4;

%Low-pass
f_cut = 70;
wc = f_cut/(fs_EEG/2);
[b,a] = butter(filter_order,wc,'low');
%figure; freqz(b1,a1,1024,fs_EEG);

figure 
[hz,hp,hl]=zplane(b,a);

set(hz, markerSize=25, color='b',linewidth=2);
set(hp, markerSize=25, color='r',linewidth=2);
set(hl,linewidth=2)
title('Low Pass Filter (70 Hz)','FontSize',16)


EEG = filtfilt(b,a,EEG);

%%
%High-pass
f_cut = 0.5;
wc = f_cut/(fs_EEG/2);
[b,a] = butter(filter_order,wc,'high');
%figure; freqz(b1,a1,1024,fs_EEG);
figure 
[hz,hp,hl]=zplane(b,a);

set(hz, markerSize=25, color='b',linewidth=2);
set(hp, markerSize=25, color='r',linewidth=2);
set(hl,linewidth=2)
title('High Pass Filter (0.5 Hz)','FontSize',16)


EEG = filtfilt(b,a,EEG);

%% Plot 
t= (1:10*fs_EEG)/fs_EEG;
figure 
plot(t,EEG(3600*fs_EEG:3610*fs_EEG-1), 'Color', '#123356');

xlabel('Time [s]');
ylabel('Amplitude [mV]')

%%
%plotting power spectrum
noverlap = 0.5*N;
[Pxx, freq1] = pwelch(EEG, hanning(length(EEG)), noverlap, [], fs_EEG);

figure
plot(freq1, Pxx)
xlim([0.5 70])

%% 

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

%% To see the result
figure
plot((1:30*fs_EEG)/fs_EEG, ep{689})
xlim([10,20])
xlabel('[s]')
ylabel('[mV]')

figure
plot((1:30*fs_EEG)/fs_EEG, ep{46})
xlim([10,20])
xlabel('[s]')
ylabel('[mV]')
%%
%Plot of the variances for each frequency band
lim = [0 9];
figure
t = (1:lep:lep*floor(length(EEG)/ep_samples))/3600;

plot(t,var.alpha),ylim([-10 50]), xlim(lim);
hold on
plot(t,var.theta),ylim([-10 50]), xlim(lim);
hold on 
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
                domTab{i,3} = 5; %awake
            end
        end

        if b>=3 && a >= c
            for j = i:i+5
                domTab{i,3} = 4; %REM
            end
        end

        if b >= 3 && d >= e
            for j = i:i+5
                domTab{i,3} = 2; %NREM 2 sleep
            end
        end

        if c >= 3 && d >= 5
            for j = i:i+5
                domTab{i,3} = 1; %NREM 3 deep sleep
            end
        end
        

      
    end

    if isempty(domTab{i, 3}) %NREM 1 light sleep
        domTab{i,3} = 3;
    end

end

%% Per fare l'hypnogram più leggibile 
app=cell2mat(domTab(:,3));
EpTogether = 10;
for i = 1:EpTogether:length(app)-EpTogether
    count=zeros(5,1);
    for k=1:5
        for j=i:i+EpTogether-1
            if app(j)==k
                count(k)=count(k)+1;
            end
        end     
    end
    [pp,p]=max(count);
    for j=i:i+EpTogether-1
        app(j,1)=p;
    end
end

%%
figure
for h = 1:1060
    hypno(h) = domTab{h,3};
end

plot(1:1060, hypno, 'LineWidth',0.75,'Color','red')


for h = 1:1060
    hypno(h) = app(h);
end

hold on
plot(1:1060, hypno,'LineWidth',3,'Color','Blue')
hold on
yline(2.5,'k--','Slow Wave Sleep','LineWidth',1)
title('Numero di Epoche insieme: ',EpTogether)
xlabel('Epochs','FontSize',15)
ylabel('Sleep Stages', 'FontSize',15)
set(gca,'YTick',0:5);
yticklabels({'';'NREM3';'NREM2';'NREM1';'REM';'AWAKE'})
ylim([0,5]);


%%
domTab = cell2table(domTab, 'VariableNames', ["Less Dominant Band", "Dominant band", "Stage"]);
save("dom_waves.mat", "domTab");
