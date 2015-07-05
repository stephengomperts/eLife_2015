%Compute VTA unit coordination with hippocampal theta

% This code computes the hippocampal theta phase associated with each VTA spike
% and the mean preferred phase and circular concentration of each VTA unit at hippocampal theta. 

%prerequisite data:
% BlackmanFilter4_12HzCoefficients :   vector of 500 coefficients generated using the filter design tool (Fs 1000 Hz)
% lfp :                                struct array containing fields of time (sec) and data (hippocamal local field potential (microvolts)) downsampled to 1KHz.  
% startime :                           start time of recording session (sec)
% stoptime :                           stop time of recording session (sec)
% VTAunits :                           struct array containing the spike times of each VTA unit recorded within a session (sec)
% time :                               (sec)
% runvelocity :                        (cm/sec)

%step 1:  from the hippocampal local field potential, compute the instantaneous theta phase

thetabin.data = filtfilt(BlackmanFilter4_12HzCoefficients,1,lfp.data);
thetabin.timestamp = lfp.timestamp;
hil = hilbert(thetabin.data);
thetabin.phase = (angle(hil));

a=[];
b=[];
spikes = [];
VTAspeed = [];
temp = [];
for i=1:size(VTAunits,2) %keep VTA unit spikes present during run behavior
    a = find(VTAunits(:,i).spikes > startime, 1, 'first');
    b = find(VTAunits(:,i).spikes < stoptime, 1, 'last');
    spikes =  VTAunits(:,i).spikes(a:b);
    VTAspeed = abs(interp1(time, runvelocity, spikes, 'linear')) ;
    temp = VTAspeed >= 10;
    VTAunits(:,i).RUNthresh = spikes(temp);
end

%step 2: for each VTA unit, compute the rayleigh statistic, mean phase and
%circular concentration coefficient at hippocampal theta

prayleigh = [];
mu = [];
kappa = [];
for i=1:size(VTAunits,2) 
    spikephase(i).VTA = interp1(thetabin.timestamp, thetabin.phase(1:numel(thetabin.timestamp)), VTAunits(i).RUNthresh, 'linear');
    prayleigh(i).VTA = rayleigh(spikephase(i).VTA);
    [mu(i).VTA,kappa(i).VTA]=vonmises_fit(spikephase(i).VTA);
end

