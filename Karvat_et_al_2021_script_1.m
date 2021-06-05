% This script shows the calculations performed to obtain Figure 2 with raw
% data from one session
clearvars -except LFP FOM;
% to calculate power we use the fieldtrip toolbox, which has to be added to the path. 
% Can be downloaded from https://www.fieldtriptoolbox.org/
addpath('D:\fieldtrip'); ft_defaults; 

% Load saved data
gdir = 'D:\THE_GitHub\';    % the folder with saved data from GitHub
load ([gdir 'FR']);         % population firing rate
load ([gdir 'ratP']);       % structure with behavioral results and timestamps per trial
if ~exist('LFP')==1         % raw LFP data from one shank (16 channels)
    for li = 1:3; load ([gdir 'lfp' int2str(li)]); end
    LFP = [lfp1,lfp2,lfp3]; % the original LFP file was split into 3 parts to allow uploading to GitHub
    clear lfp1 lfp2 lfp3 li
end

% parameters
fs          = 976.5625;                                 % sampling frequency in Hz
foi         = 1:120;                                    % frequencies of interest, in Hz
band_lims   = {[3 10],[15 30],[45 90],[95 120]};        % The borders of bands, in Hz
band_names  = {'Low','Beta','Low gamma','High gamma'};
vib_ON      = [ratP([ratP.res_num]==1).piezo_ON];       % the time (sample) in which the (real) piezo was activated in Hit trials
ctrl_ON     = [ratP([ratP.res_num]>2).PiezoCtrl_ON];    % the time (sample) in which the control piezo was activated 

%% ERP
pre         = 0.05;                                         % time is sec before stimulus onset
post        = 0.15;                                         % time is sec after stimulus onset
inds        = floor(-pre*fs) : ceil(post*fs);
t           = inds/fs;                                      % time vector
siz         = [size(LFP,1),size(inds,2),size(vib_ON,2)];    %dimord: channel_time_trial
LFP_trialed = zeros(siz);
for tr = 1:length(vib_ON)                                   % divide into trials
    LFP_trialed(:,:,tr) = LFP(:,inds+ vib_ON(tr));
end
ERP = mean(LFP_trialed,3);                                  % the ERP is the mean of LFP over trials
figure(111); clf; % plot
imagesc(t, 1:siz(1), ERP);
colormap (jet); colorbar
xlabel ('Time (sec)'); ylabel('Channel #'); title('ERP');

%% CSD
a = LFP_trialed(1:end-2,:,:);   % channels above
b = LFP_trialed(3:end,:,:);     % channels below
c = LFP_trialed(2:end-1,:,:);   % target channel
CSD = (a + b - 2*c);            % 2nd spatial derivative estimate
clear a b c
meanCSD = mean(CSD,3);          % mean over trials
figure(112); clf; % plot
imagesc(t, 1:siz(1), meanCSD);
colormap (jet); colorbar
xlabel ('Time (sec)'); ylabel('Channel #'); title('CSD');

%% Population firing rate 
toi     = [-0.5 0.5];                           % time (sec) of interest to plot
inds    = floor(toi(1)*fs) : ceil(toi(2)*fs);
t       = inds/fs;                              % time vector
FR_Hit  = FR(inds+vib_ON');                     % divide into Hit trials
FR_Ctrl = FR(inds+ctrl_ON');                    % divide into control trials

figure(121); clf; hold on % plot
plot(t, mean(FR_Ctrl), 'color',ones(1,3)*0.5); 
plot(t, mean(FR_Hit), 'b'); 
title('Population firing-rate');
xlim(toi);
xlabel('Time (sec)');
ylabel ('FR (z-score)');

%% Induced power 
% to allow working in most computers, the TFR is calculated channel-by-channel
if ~(exist('FOM')==1)
FOM = zeros(size(LFP,1), length(foi), size(LFP,2)); % matrix for fraction of the median of power in all channels
for chi = 1:size(LFP,1)                             % on a regular PC this loop takes ~1 minute per channel
    % make the LFP a fieldtrip structure
    cfg          = [];
    cfg.raw_data = LFP(chi,:);
    cfg.SF       = fs;
    cfg.demean   = 'no';
    cfg.detrend  = 'no';
    cfg.resample = 0;
    LFP_raw = data2structGH(cfg);
    % make the TFR
    cfg             = [];
    cfg.method      = 'wavelet';    %time-frequency analysis using the 'wavelet method' based on Morlet wavelets
    cfg.output      = 'pow';        % the output wll be power
    cfg.foi         = foi;
    cfg.toi         = 'all';        % run the TFR analyis over the entire session
    cfg.width       = 7;            % the Morlet wavelet will have width of 7 cycles for each frequencies
    cfg.pad         = 'nextpow2';   % makes analysis faster
    TFR = ft_freqanalysis(cfg, LFP_raw);
    pwr = squeeze(TFR.powspctrm);   % take only the poewr out of the fieldtrip structure
    
    cfg.art_thresh  = 1000;         % threshold (in uV) for artifact removal
    cfg.art_removal = 500;          % duration (in samples) before and after an artifact to be removed (put NaN).
    pwrClean = gross_artifact_removal_GH(cfg, LFP(chi,:), pwr); % artifact removal. Crucial for percentile calculation
    pwr_ratio = pwrClean./nanmedian(pwrClean,2); % Normalize to the median (FOM)
    FOM(chi,:,:) = pwr_ratio;
end
clear pwr_ratio pwrClean
end

toi =       [-0.5 0.5];                                             % time (sec) of interest to plot
inds =      floor(toi(1)*fs) : ceil(toi(2)*fs);
t =         inds/fs;                                                % time vector
siz =       [size(LFP,1),length(foi),size(inds,2),size(vib_ON,2)];  %dimord: channel_freq_time_trial
FOM_trialed = zeros(siz);
for tr = 1:length(vib_ON) % divide into trials
    FOM_trialed(:,:,:,tr) = FOM(:,:,inds+ vib_ON(tr));
end

figure (131); clf; % plot
for bi = 1:4 % per band
    curr = FOM_trialed(:,band_lims{bi}(1):band_lims{bi}(2),:,:); % take the current band
    curr = squeeze(nanmean(curr,2)); % take mean power of the band
    curr = squeeze(nanmean(curr,3)); % average over trials    
    subplot(4,1,bi);
    imagesc(t, 1:size(curr,1), curr);
    colorbar
    title(band_names{bi}); 
    if bi==4; xlabel('Time (sec)'); end
    ylabel ('Electrode #');
end
colormap jet

%% Burst analysis
% For example, beta bursts on channel 5
pctl         = 90; % the power percentile above which bursts should be

% make the LFP a fieldtrip structure
cfg          = [];
cfg.raw_data = LFP(5,:);
cfg.SF       = fs;
cfg.demean   = 'no';
cfg.detrend  = 'no';
cfg.resample = 0;
LFP_raw = data2structGH(cfg);
% make the TFR
cfg             = [];
cfg.method      = 'wavelet';                             % time-frequency analysis using the 'wavelet method' based on Morlet wavelets
cfg.output      = 'pow';                                 % the output wll be power
cfg.foi         = band_lims{2}(1)-1 : band_lims{2}(2)+1; % take beta, plus 1 Hz at edges to allow detection of peaks 
if cfg.foi(1)==0; cfg.foi(1) = []; end    
cfg.toi         = 'all';                                 % run the TFR analyis over the entire session
cfg.width       = 7;                                     % the Morlet wavelet will have width of 7 cycles for each frequencies
cfg.pad         = 'nextpow2';                            % makes analysis faster
TFR = ft_freqanalysis(cfg, LFP_raw);
pwr = squeeze(TFR.powspctrm);

cfg.art_thresh  = 1000;                                  % threshold (in uV) for artifact removal
cfg.art_removal = 500;                                   % duration (in samples) before and after an artifact to be removed (put NaN).
pwrClean = gross_artifact_removal_GH(cfg, LFP(5,:), pwr);% artifact removal. Crucial for percentile calculation
pwr_masked = pwrClean./prctile(pwrClean,pctl,2);         % normalize the clean power to the (90th) percentile
pwr_masked(pwr_masked < 1)=0;                            % put a mask over all values less than the (90th) percentile
pwr_masked(isnan(pwr_masked))=0;                         % artifacts had NaN, which have to be converted to 0
BW = imregionalmax(pwr_masked);                          % detect peaks in the 2D plane
BW ([1,end],:) = [];                                     % one Hz above and one below were added just for peak reasons
[i,j] = find(BW);                                        % get the coordinates of peaks (bursts)

burstDists  = zeros(length(j),4);   %1- freq (Hz), 2- timing (sample), 3- begin, 4- end
bursVec     = zeros(1,size(LFP,2)); % a vector with "1" when a burst is "on", and "0" when "off"
badPeaks    = [];                   % to remove bursts with an error in calculation
t_trials    = tic;                  % for progress report
prev        = fprintf(' ');
for x = 1:length(j)
    if ~mod(x,500); prev = dispRMVprevGH (['Finished ' num2str(x) ' bursts in ' num2str(toc(t_trials)) ' sec'],prev); end
    try
        t_before = find(~pwr_masked(i(x),1:j(x)),1,'last')+1; if isempty(t_before); t_before=1; end
        t_after  = find(~pwr_masked(i(x),j(x):end),1,'first') + j(x) - 2; if isempty(t_after); t_after = size(pwr_masked,2); end
        burstDists(x,1) = i(x)+cfg.foi(1);      % peak frequency
        burstDists(x,2) = j(x);                 % peak timing
        burstDists(x,3) = t_before;             % samples before
        burstDists(x,4) = t_after;              % samples after
        bursVec(t_before:t_after) = 1;          % update the burst vector
    catch
        badPeaks = [badPeaks, x];               % log if there was an error in burst detection
    end
end
burstDists(badPeaks,:) = [];

toi     = [-0.5 0.5];                           % time (sec) of interest to plot
inds    = floor(toi(1)*fs) : ceil(toi(2)*fs);
t       = inds/fs;                              % time vector
BO_Hit  = bursVec(inds+vib_ON');
mpwr    = pwrClean(2:end-1,:);                  % to compare to power, take the beta frequencies
mpwr    = mpwr./nanmedian(mpwr,2);              % normalize to the median (FOM)
mpwr    = nanmean(mpwr);                        % take the mean of beta frequencies
pwr_Hit = mpwr(inds+vib_ON');                   % divide into trials

figure (141); clf; hold on % plot
yyaxis left
plot(t, mean(BO_Hit));
ylabel ('Burst occupancy');
yyaxis right
plot(t, mean(pwr_Hit));
ylabel ('Beta power (FOM)');
title('Beta burst occupancy and power');
xlim(toi);
xlabel('Time (sec)');
%% Beta burst localization
pre     = 0.2;                                 % time is sec before stimulus onset
post    = 0.2;                                 % time is sec after stimulus onset
inds    = floor(-pre*fs) : ceil(post*fs);
t       = inds/fs;                             % time vector
siz     = [size(FOM,1), size(FOM,2),  length(t)]; % dimord: chan_freq_time
binds = burstDists(:,2)+inds;
BstLoc = zeros(siz);
for chi = 1:siz(1)
    for fi = 1:siz(2)
        curr = squeeze(FOM(chi,fi,:));
        curr = nanmean(curr(binds));
        BstLoc(chi,fi,:) = curr;
    end
    disp (['Finished channel ' num2str(chi)]);
end
%%
clim = [min(BstLoc(:)), max(BstLoc(:))];
max_f = 50;
figure(151); clf;
for chi = 1:siz(1)
    subplot(siz(1),1,chi);
    imagesc(t,1:max_f,squeeze(BstLoc(chi,1:max_f,:)),clim);
    axis xy
    text(0.05, 0.5, num2str(chi), 'fontweight','bold','units','normalized','color','w');
    if chi==siz(1); xlabel('Time (sec)'); end
    if chi==1; title('beta burst localization'); end
end
colormap jet