% clear all;

load('./participant2.mat');



[b,a]=butter(4,[5,10]/(1000/2),'bandpass');
% fullecg=-data(3*fs:1850*fs,1);
fullecg=-data(:,1);
% fullecg=-data(4*fs:1847*fs,1);
yECG=filtfilt(b,a,fullecg);
t = (0:length(yECG(:,1))-1)/fs;
plot(t,yECG);

%%

participant2_times2;

sections = {'rest1', 'mistc', 'mist1', 'mist2', 'mist3', ...
            'rest2', 'stroop1', 'stroop2', 'restandmistc', 'mist123'};

% sections = {'restandmistc', 'mist123'};

all = {'meanRR', 'stdRR',    ...
       'meanHR', 'stdHR',    ...
       'RMSSD', 'SDSD',      ...
       'NN50', 'PNN50',      ...
       'VLF', 'HF', 'LF',    ...
       'PTOT', 'LFN', 'HFN', ...
       'LFHF',               ...
       'SampEn1', 'SampEn2', ...
       'meanIBI', 'stdRR', ...
       'meanResp', 'stdResp'
       };

%%
% 
relax_seg = fs*t_before_relax:fs*(t_before_relax+180);
mistc_seg = fs*t_before_mist_c:fs*(t_before_mist_c+180);
mist1_seg = fs*t_before_mist_1:fs*(t_before_mist_1+180);
mist2_seg = fs*t_before_mist_2:fs*(t_before_mist_2+180);
mist3_seg = fs*t_before_mist_3:fs*(t_before_mist_3+180);
relax2_seg = fs*t_before_relax_2:fs*(t_before_relax_2+180);
stroop1_seg = fs*t_before_stroop_1:fs*(t_before_stroop_1+180);
stroop2_seg = fs*t_before_stroop_2:fs*(t_before_stroop_2+180);
restandmistc_seg = relax_seg(1):mistc_seg(end);
mist123_seg = mist1_seg(1):mist3_seg(end);



raw_ecg = {data(relax_seg,1);   ...
           data(mistc_seg,1);   ...
           data(mist1_seg,1);   ...
           data(mist2_seg,1);   ...
           data(mist3_seg,1);   ...
           data(relax2_seg,1);   ...
           data(stroop1_seg,1); ...
           data(stroop2_seg,1); ...
           data(restandmistc_seg,1);   ...
           data(mist123_seg,1);   ...
           };
       
raw_resp ={data(relax_seg,2);   ...
           data(mistc_seg,2);   ...
           data(mist1_seg,2);   ...
           data(mist2_seg,2);   ...
           data(mist3_seg,2);   ...
           data(relax2_seg,2);   ...
           data(stroop1_seg,2); ...
           data(stroop2_seg,2); ...
           data(restandmistc_seg,2);   ...
           data(mist123_seg,2);  ...
           };
       
% raw_ecg = {data(restandmistc,1);   ...
%            data(mist123,1);   ...
%            };
%        
% raw_resp ={data(restandmistc,2);   ...
%            data(mist123,2);   ...
%            };

% for s = 8:length(sections)
for s = 1:length(sections)
    pause(3);
    fprintf('Next: %s.\n', sections{s});
    input('ready for next?');
    close all;
    keep raw_ecg raw_resp sections data fs s all all_stats altpow
    [xRRI,fsRRI,pECG,tpECG]=ECG_to_RRI(cell2mat(raw_ecg(s)), fs);
    % xRRI:    interpolated RRI
    % fsRRI:   sample rate for interpolation
    % pECG:    ECG peaks (peakECG)
    % tpECG:   time of ECG peaks (time-peakECG)
    % dtpECG:  difference of ECG peak times (difference-time-peaksECG)
    % fdtpECG: filtered difference of ECG peak times
    %          (filtered-difference-time-peaksECG).
    
    % Time Domain Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meanRR = mean(xRRI);         % (ms)
    stdRR = sqrt(var(xRRI));     % (ms)
    meanHR = 60/meanRR;          % (mins^-1)
    stdHR = sqrt(var(60./xRRI)); % (mins^-1)

    dtpECG = diff(tpECG);
    % Remove outliers in RR timests due to rem. of outliers in RR series.
    tol = 2; % Only accept RR intervals that 
             % are less than 2x the mean.
    fdtpECG = dtpECG( dtpECG < tol*mean(dtpECG)); % filtered peaks of ECG

    RMSSD = sqrt(mean(diff(fdtpECG).^2)); % (ms)
    SDSD = std(fdtpECG);                  % (ms)

    NN50 = sum(diff(fdtpECG) > 0.050);    % (quantity)
    PNN50 = NN50/length(dtpECG);         % (%)

    % NN20, suggested by Mietus et. al. (2002)
    NN20 = sum(diff(fdtpECG) > 0.020);    % (quantity)
    PNN20 = NN20/length(dtpECG);         % (%)
    
    % Frequency Domain Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [pxx_ls, f_ls] = plomb(pECG,tpECG);

%     [d,p] = aryule(xRRI-mean(xRRI),16);

    [ppxx, pf] = pyulear(xRRI(1:end-50)-mean(xRRI(1:end-50)),16, 2^15, fsRRI);

    pVLF_idx = (pf>=0) & (pf <= 0.04);
    pLF_idx = (pf>=0.04) & (pf <= 0.15);
    pHF_idx = (pf>=0.15) & (pf <= 0.40);

    VLF = sum(ppxx(pVLF_idx));
    LF = sum(ppxx(pLF_idx));
    HF = sum(ppxx(pHF_idx));

    
%     VLF_idx = (pf>=0) & (f_ls <= 0.04);
%     LF_idx = (f_ls>=0.04) & (f_ls <= 0.15);
%     HF_idx = (f_ls>=0.15) & (f_ls <= 0.40);
%     
%     VLF = sum(pxx_ls(VLF_idx));
%     LF = sum(pxx_ls(LF_idx));
%     HF = sum(pxx_ls(HF_idx));
    PTOT = sum(ppxx(pVLF_idx|pLF_idx|pHF_idx));
    
%     PTOT = sum(pxx_ls(VLF_idx|LF_idx|HF_idx));
    
    LFN = 100* (LF/(PTOT-VLF));
    HFN = 100* (HF/(PTOT-VLF));
    
    LFHF = LF/HF;
    
    % Nonlinear Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N=length(xRRI);
    sd = std(xRRI);
    r = 0.15;
    
    norm_xRRI = zscore(xRRI');
    
    tau=1;
    clear y
    for j = 1:N/tau
        y(j) = (1/tau) * sum(xRRI((j-1)*tau+1:j*tau));
    end
    e = sampen(y, 2, r, 0, 0, 0);
    SampEn1 = e(2);

    tau=2;
    clear y
    for j = 1:N/tau
        y(j) = (1/tau) * sum(xRRI((j-1)*tau+1:j*tau));
    end
    e = sampen(y, 2, r, 0, 0, 0);
    SampEn2 = e(2);
    
    % Respiration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IBI = cjw_resp(cell2mat(raw_resp(s)),fs);
    
    meanIBI = mean(IBI);           % (ms)
%     stdIBI = sqrt(var(IBI));        % (ms)
    meanResp = 60/IBI;          % (mins^-1)
%     stdResp = sqrt(var(60./IBI));  % (mins^-1)
    
    [pxx, f] = periodogram(xRRI,hamming(length(xRRI)), 2.^(2+nextpow2(length(xRRI))), fsRRI);
    VLFpow = bandpower(pxx,f,[0 0.04],'psd');
    LFpow = bandpower(pxx,f,[0.04 0.15],'psd');
    HFpow = bandpower(pxx,f,[0.15 0.4],'psd');
    
    altpow(s,1) = VLFpow;
    altpow(s,2) = LFpow;
    altpow(s,3) = HFpow;
    altpow(s,4) = LFpow/HFpow;
    
    for iii = 1:length(all)
         all_stats(s,iii) = eval(all{iii});
    end
end

%%

% [H,w]=freqz(sqrt(p),d);
% plot(w, H);
% xlim([0 0.5]);

[pxx_w, f_w] = pwelch(xRRI-mean(xRRI),hamming(round(length(xRRI)./2)),round(length(xRRI)./4), 2.^(nextpow2(length(xRRI))), fsRRI);
% [d,p] = aryule(xRRI-mean(xRRI),16);
% [pxx_w, f_w] = pyulear(xRRI,16, 1024, fsRRI);

plot(f_w, pxx_w);
xlim([0 0.5]);
ylim([0 0.03]);

    pVLF_idx = (f_w>=0) & (f_w <= 0.04);
    pLF_idx = (f_w>=0.04) & (f_w <= 0.15);
    pHF_idx = (f_w>=0.15) & (f_w <= 0.40);

    VLF = sum(pxx_w(pVLF_idx));
    LF = sum(pxx_w(pLF_idx));
    HF = sum(pxx_w(pHF_idx));

% VLFpow = bandpower(pxx_w,f_w,[0 0.04],'psd')
% LFpow = bandpower(pxx_w,f_w,[0.04 0.15],'psd')
% HFpow = bandpower(pxx_w,f_w,[0.15 0.4],'psd')
