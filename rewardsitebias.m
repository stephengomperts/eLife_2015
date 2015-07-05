%Computing reward site bias

% This code computes a reward site bias for each replay time bin
% and for each VTA spike-associated replay time bin.

%prerequisite data:
% SPWReventstartime :      vector of SPW-R event start times (sec)
% replayeventstartime :    vector containing the subset of SPWReventstartime that are replay events and that happen at a specified location (for example, at Force reward sites)
% bins10 :                 vector of the 10cm spatial bins that span the track
% VTAunits :               struct array containing the spike times of each VTA unit recorded within a session (sec)
% BINSIZE :                replay event time bin size (25 msec)
% pdfall :                 array of spatially decoded, concatenated candidate replay events, comprising a spatial pdf (10cm bins) x time (25msec bins) 
% candidatebinlength :     vector of SPW-R time bin lengths
% Rsites :                 spatial bins associated with reward sites
% armsites :               spatial bins not associated with reward sites

RsiteBIAS=[];
VunitRsiteBIAS = [];

for i=1:length(SPWReventstartime)-1;
    estA=[];
    
    %step 1: across SPWR events, select replay events occurring at the specified location, to
    % derive pdf's of position and direction
    
    if intersect(SPWReventstartime(i), replayeventstartime)>0
        xi = sum(candidatebinlength(1:i))+1;
        yi = sum(candidatebinlength(1:i+1));
        replaybins = sum(candidatebinlength(1:i))+1:sum(candidatebinlength(1:i+1));
        estA = pdfall(bins10,xi:yi); %single replay event spatial pdf x time
        
        %step 2: across each time bin of the replay event, derive reward site bias
        
        for k = 1:length(replaybins)
            Rcon = [];
            armcon = [];
            RconAVG = [];
            armconAVG = [];
            
            Rcon=estA(Rsites,k);
            armcon=estA(armsites,k);
            RconAVG = (mean(estA(Rsites,k)));
            armconAVG = (mean(estA(armsites,k)));
            RsiteBIAS = [RsiteBIAS (RconAVG-armconAVG)];
        end
        
        %step 3: across replay bins and VTA units, acquire replay bin spatial content associated with VTA unit spikes
        
        for k = 1:length(replaybins)
            VunitRconAVG = [];
            VunitarmconAVG = [];
            Vunitpos = zeros(size(bins10,2),size(VTAunits,2));
            delay =  0.084;
            x10 = binT10(SPWReventstartime(i))+(k-1)*BINSIZE + delay;
            y10 = binT10(SPWReventstartime(i))+(k*BINSIZE) + delay;
            
            %step 4: across VTA units, weight replay bin spatial content by the number of associated VTA unit spikes
            
            for z=1:size(Vunitpos,2) 
                weight =(length(find(VTAunits(z).spikes>=x10,1,'first'):find(VTAunits(z).spikes<y10,1,'last')));
                if weight>0
                    Vunitpos(:,z) = Vunitpos(:,z) + ...
                        estA(:,k)*weight;
                end
            end
            
            VunitRconAVG = (mean(Vunitpos(Rsites,:)));
            VunitarmconAVG = (mean(Vunitpos(armsites,:)));
            
            %account for 0 VTA spikes at a replay bin
            VunitRconAVG(VunitRconAVG==0)=nan;
            VunitarmconAVG(VunitarmconAVG==0)=nan;
            
            VunitRsiteBIAS = [VunitRsiteBIAS (VunitRconAVG-VunitarmconAVG)'];
        end
    end
end

%step 5: transform Rsitebias and VunitRsitebias into binary output

RsiteBIAS(RsiteBIAS==0) = nan;
RsiteBIAS(RsiteBIAS>0) = 1;
RsiteBIAS(RsiteBIAS<0) = 0;

VunitRsiteBIAS(VunitRsiteBIAS>0) = 1;
VunitRsiteBIAS(VunitRsiteBIAS<0) = 0;



