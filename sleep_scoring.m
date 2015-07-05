%Sleep scoring

% This code identifies sleep segments and sleep-associated SPW-R events, 

%prerequisite data:
% HCmultiunitspiketimes :                  vector of hippocampal multiunit spike times (sec)
% BlackmanFilter4_12HzCoefficients :       vector of coefficients generated using the filter design tool (Fs 1000 Hz)
% BlackmanFilter_under4Hz_Coefficients :   vector of coefficients generated using the filter design tool (Fs 1000Hz)
% startime, stoptime :                     start and stop time of sleep session (sec)
% lfp :                                    struct array containing fields of time (sec) and data (hippocampal local field potential (microvolts)) downsampled to 1KHz.  
% emg  :                                   struct array containing fields of time (sec) and data (emg (microvolts)) downsampled to 1KHz. 

%step 1: compute mean theta power/delta power ratio and mean EMG power in 2 second bins across the sleep session

binT = (startime:0.010:stoptime); 
PCLmu_bin = histc(HCmultiunitspiketimes, binT); %bin all spikes, from startime to stoptime
smPCLmu_bin = smoothn(PCLmu_bin, 1); %smooth with a running gaussian window of width 1 sample (10ms)
theta = filtfilt(BlackmanFilter4_12HzCoefficients,1,lfp.data); %Fs 1000
delta = filtfilt(BlackmanFilter_under4Hz_Coefficients,1,lfp.data); %Fs 1000
EMG = emg.data;

thetaint = [];
deltaint = [];
tdratioint = [];
EMGint = [];
time = [];
bin = 2;
j=0;
for i = bin:bin:fix(stoptime-startime) 
    j=j+1;    
    t1 = find(lfp.time >startime-bin+i,1,'first');
    t2 = find(lfp.time >=startime+i,1,'first');
    thetaint(j) = mean(theta(t1:t2).^2) ;
    deltaint(j) = mean(delta(t1:t2).^2) ;
    EMGint(j) = mean(EMG(t1:t2).^2) ;
end
tdratioint = thetaint/deltaint;

%step 2: compute sleep-associated candidate events as peaks in the smoothed hippocampal multiunit

segs_SPWRevent = detect_mountains(binT, smPCLmu_bin, ...
    'threshold', [mean(smPCLmu_bin),mean(smPCLmu_bin)+4*std(smPCLmu_bin)], ...
    'span_interval', 0.020, ...
    'width_lim', 0.020);

%step 3: find peaks in the theta/delta ratio

segsTD = detect_mountains(time, tdratioint, ...
    'threshold', nanmean(tdratioint), ...
    'span_interval', bin, ...
    'width_lim', bin);

%step 4: find peaks in the EMG power

segsEMG = detect_mountains(time, EMGint, ...
    'threshold', nanmean(EMGint), ...
    'span_interval', bin, ...
    'width_lim', bin);

%step 5: select slow wave sleep segments (array of slow wave sleep segment start and stop times)

lowTD = seg_excl([startime stoptime],segsTD); % exclude segments with high theta/delta ratio

segsSWS = seg_excl(lowTD,segsEMG); % exclude segments with high EMG power (associated with movement)

%step 6: select slow wave sleep-associated SPW-R events (array of SPW-R event start times and stop times)

SPWRevent_SWS = seg_excl(segs_SPWRevent,segsTD); %exclude segments with high theta/delta segmemts

SPWRevent_SWS = seg_excl(SPWRevent_SWS,segsEMG); %exclude segments with high EMG power




















