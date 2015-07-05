%Compute reward site bias for subtypes of replay

% This code computes a reward site bias for each replay time bin and for each 
% VTA spike-associated replay time bin, restricting the analysis to either 
% (a) outbound or inbound replay time bins ('direction' in track coordinates), 
% and to replay events with either (b) outgoing (a positive replay slope away
% from) or incoming (a negative replay slope towards) the rat's location.
% Forward replay occurs with (a) the intersection of positive replay slope and outbound
% direction, and with (b) the intersection of negative replay slope and inbound direction. 
% Reverse replay occurs with (a) the intersection of positive replay slope and inbound
% direction, and with (b) the intersection of negative replay slope and outbound direction.
% Centrifugal replay occurs with positive replay slope for both directions. 
% Centripetal replay occurs with negative replay slope for both directions.

%prerequisite data:
% SPWReventstartime :      vector of SPW-R event start times (sec)
% replayeventstartime :    vector containing the subset of SPWReventstartime 
%                          that are replay events and that happen at a specified 
%                          location (for example, at Force reward sites on the SWM task)
% bins10 :                 vector of the 10cm spatial bins that span the track
% VTAunits :               struct array containing the spike times of each VTA
%                          unit recorded within a session (sec)
% BINSIZE :                replay event time bin size (25 msec)
% pdfall :                 array of spatially decoded, concatenated candidate 
%                          replay events, comprising a spatial pdf (10cm bins) x time (25msec bins)
% pdfdirection :           array of direction decoded, concatenated candidate 
%                          replay events, comprising a direction pdf (outbound/inbound) x time (25ms bins)
% candidatebinlength :     vector of SPW-R time bin lengths
% Rsites :                 vector of spatial bins associated with reward sites
% armsites :               vector of spatial bins not associated with reward sites
% replayslope              vector of replay event slopes (m/sec), computed using the
%                          radon analysis. See Davidson and Kloosterman, 2009
%                          and Kloosterman, 2012 for instructions. 
% delay :                  latency between the hippocampal SPW-R and the VTA (0.084 seconds)
% DIRthresh :              threshold for replay direction preference (0.6)        
% positionatreplay :       vector of the animal's location (cm) at each replay event
% select_tracktype :       user specified string to distinguish between a recording on the SWM task or on the linear track ('SWMtrack' or 'lineartrack')
% select_direction :       user specified string to acquire either outbound ('POS') or inbound ('NEG') direction-associated reward site bias
% select_replayslope :     user specified string to acquire either outgoing replay slope ('outgoing') or incoming replay slope ('incoming')



%step 1: assign select_direction and select_replayslope, then proceed as follows.

RsiteBIASPOS = [];
RsiteBIASNEG = [];
if strcmp(select_replayslope,'outgoing') % select outgoing replay events
    for i=1:length(SPWReventstartime)-1;  % across all SPWR events
        estA = [];
        estD = [];
        
        % step 2: select replay events occurring at the specified location (Force
        % reward locations), to derive pdf's of position and direction
        
        if intersect(SPWReventstartime(i), replayeventstartime) 
            Vunitpos_nothresh = zeros(size(bins10,2),size(VTAunits,2));
            xi = sum(candidatebinlength(1:i))+1;
            yi = sum(candidatebinlength(1:i+1));
            replaybins = sum(candidatebinlength(1:i))+1:sum(candidatebinlength(1:i+1));
            
            %single replay event spatial pdf x time 
            estA = pdfall10(bins10,xi:yi);  
            
            %single replay event direction pdf x time 
            estD = pdfdirection10(:,xi:yi); 
            
            %inbound direction estimate
            inbound = estD(1,:);  
            
            %outbound direction estimate
            outbound = estD(2,:); 
            
            pdfdirection_index = (outbound-inbound)./(outbound+inbound);
            
            estAPOS = zeros(size(estA));
            estANEG = zeros(size(estA));
            
            if ~strcmp(select_tracktype, 'linear') %consider SWM task
                
                %for outgoing replay events with rat at force reward sites 
                %(outbound direction estimate = forward replay, inbound direction estimate = reverse replay)
                
                if replayslope(i) > 0 
                    
                    %step 3: across replay bins, compute POS and NEG reward site bias in each bin
                    
                    for k = 1:length(replaybins) 
                        RconAVGPOS = [];
                        armconAVGPOS = [];
                        
                        RconAVGNEG = [];
                        armconAVGNEG = [];
                        
                        estAPOS(:,k) = estA(:,k) * (pdfdirection_index(k) > DIRthresh);   
                        estANEG(:,k) = estA(:,k) * (pdfdirection_index(k) < -DIRthresh);  
                        
                        RconAVGPOS = (mean(estAPOS(Rsites,k)));
                        armconAVGPOS = (mean(estAPOS(armsites,k)));
                        RsiteBIASPOS = [RsiteBIASPOS (RconAVGPOS-armconAVGPOS)]; 
                        
                        RconAVGNEG = (mean(estANEG(Rsites,k)));
                        armconAVGNEG = (mean(estANEG(armsites,k)));
                        RsiteBIASNEG = [RsiteBIASNEG (RconAVGNEG-armconAVGNEG)]; 
                    end
                    
                    %step 4: across replay bins, acquire replay bin spatial content associated with VTA unit spikes
                    
                    for k = 1:length(replaybins) 
                        Vunitpos_nothreshPOS = zeros(size(bins10,2),size(VTAunits,2)); 
                        Vunitpos_nothreshNEG = zeros(size(bins10,2),size(VTAunits,2)); 
                        
                        VunitRconAVGPOS = [];
                        VunitarmconAVGPOS = [];
                        
                        VunitRconAVGNEG = [];
                        VunitarmconAVGNEG = [];
                        
                        x10 = binT10(LEstartime(i))+(k-1)*BINSIZE + delay;
                        y10 = binT10(LEstartime(i))+(k*BINSIZE) + delay;
                        [maxest,Imaxest] = max(estA(:,k));
                        
                        %across VTA units
                        for z=1:size(VTAunits,2) 
                            
                            %weight replay bin spatial content by the number of associated VTA unit spikes
                            weight =(length(find(VTAunits(z).spikes>=x10,1,'first'):find(VTAunits(z).spikes<=y10,1,'last'))); 
                            if weight>0
                                Vunitpos_nothreshPOS(:,z) = Vunitpos_nothreshPOS(:,z) + ...
                                    estAPOS(:,k)*weight;
                                
                                Vunitpos_nothreshNEG(:,z) = Vunitpos_nothreshNEG(:,z) + ...
                                    estANEG(:,k)*weight;
                            end
                        end
                        
                        VunitRconAVGPOS = (mean(Vunitpos_nothreshPOS(Rsites,:)));
                        VunitarmconAVGPOS = (mean(Vunitpos_nothreshPOS(armsites,:)));
                        
                        %account for 0 VTA spikes at a replay bin
                        VunitRconAVGPOS(VunitRconAVGPOS==0)=nan;                                                
                        VunitarmconAVGPOS(VunitarmconAVGPOS==0)=nan;                                            
                        VunitRsiteBIASPOS = [VunitRsiteBIASPOS (VunitRconAVGPOS-VunitarmconAVGPOS)'];           
                        
                        VunitRconAVGNEG = (mean(Vunitpos_nothreshNEG(Rsites,:)));
                        VunitarmconAVGNEG = (mean(Vunitpos_nothreshNEG(armsites,:)));
                        
                        %account for 0 VTA spikes at a replay bin
                        VunitRconAVGNEG(VunitRconAVGNEG==0)=nan;                                                
                        VunitarmconAVGNEG(VunitarmconAVGNEG==0)=nan;                                            
                        VunitRsiteBIASNEG = [VunitRsiteBIASNEG (VunitRconAVGNEG-VunitarmconAVGNEG)'];           
                    end
                end  
                
            elseif strcmp(select_tracktype, 'linear') %consider linear track
                
                %for outgoing replay events with rat at near (left) reward site. 
                %(inbound = reverse replay; outbound = forward replay)
                
                if replayslope(i) > 0 && positionatreplay(i)< bins(fix(length(bins)/2)) 
                    
                    %step 5: across replay bins, compute POS and NEG reward site bias in each bin
                    
                    for k = 1:length(replaybins)
                        RconAVGPOS = [];
                        armconAVGPOS = [];
                        
                        RconAVGNEG = [];
                        armconAVGNEG = [];
                        
                        estAPOS(:,k) = estA(:,k) * (pdfdirection_index(k) > DIRthresh);  %keep replay bins with strong outbound direction estimate
                        estANEG(:,k) = estA(:,k) * (pdfdirection_index(k) < -DIRthresh); %keep replay bins with strong inbound direction estimate
                        
                        RconAVGPOS = (mean(estAPOS(Rsites,k)));
                        armconAVGPOS = (mean(estAPOS(armsites,k)));
                        RsiteBIASPOS = [RsiteBIASPOS (RconAVGPOS-armconAVGPOS)]; 
                        
                        RconAVGNEG = (mean(estANEG(Rsites,k)));
                        armconAVGNEG = (mean(estANEG(armsites,k)));
                        RsiteBIASNEG = [RsiteBIASNEG (RconAVGNEG-armconAVGNEG)]; 
                    end
                    
                    %step 6: across replay bins, acquire replay bin spatial content associated with VTA unit spikes
                    
                    for k = 1:length(replaybins) 
                        Vunitpos_nothreshPOS = zeros(size(bins10,2),size(VTAunits,2)); 
                        Vunitpos_nothreshNEG = zeros(size(bins10,2),size(VTAunits,2)); 
                        VunitRconAVGPOS = [];
                        VunitarmconAVGPOS = [];
                        
                        VunitRconAVGNEG = [];
                        VunitarmconAVGNEG = [];
                        
                        x10 = binT10(LEstartime(i))+(k-1)*BINSIZE + delay;
                        y10 = binT10(LEstartime(i))+(k*BINSIZE) + delay;
                        [maxest,Imaxest] = max(estA(:,k));
                        
                        %across VTA units
                        for z=1:size(VTAunits,2)
                            weight =(length(find(VTAunits(z).spikes>=x10,1,'first'):find(VTAunits(z).spikes<=y10,1,'last'))); 
                            if weight>0
                                Vunitpos_nothreshPOS(:,z) = Vunitpos_nothreshPOS(:,z) + ...
                                    estAPOS(:,k)*weight;
                                
                                Vunitpos_nothreshNEG(:,z) = Vunitpos_nothreshNEG(:,z) + ...
                                    estANEG(:,k)*weight;
                            end
                        end
                        
                        VunitRconAVGPOS = (mean(Vunitpos_nothreshPOS(Rsites,:)));
                        VunitarmconAVGPOS = (mean(Vunitpos_nothreshPOS(armsites,:)));
                        VunitRconAVGPOS(VunitRconAVGPOS==0)=nan;
                        VunitarmconAVGPOS(VunitarmconAVGPOS==0)=nan;
                        VunitRsiteBIASPOS = [VunitRsiteBIASPOS (VunitRconAVGPOS-VunitarmconAVGPOS)'];            
                        
                        VunitRconAVGNEG = (mean(Vunitpos_nothreshNEG(Rsites,:)));
                        VunitarmconAVGNEG = (mean(Vunitpos_nothreshNEG(armsites,:)));
                        VunitRconAVGNEG(VunitRconAVGNEG==0)=nan;
                        VunitarmconAVGNEG(VunitarmconAVGNEG==0)=nan;
                        VunitRsiteBIASNEG = [VunitRsiteBIASNEG (VunitRconAVGNEG-VunitarmconAVGNEG)'];           
                    end
                    
                %consider the case of incoming decoded replay events with rat at far (right) reward site. 
                %(Outbound = reverse replay; Inbound = forward replay)
                
                elseif replayslope(i) < 0 && positionatreplay(i) > bins(fix(length(bins)/2))
                    
                    %step 7: across replay bins, compute POS and NEG reward site bias in each bin
                    
                    for k = 1:length(replaybins) 
                        RconAVGPOS = [];
                        armconAVGPOS = [];
                        
                        RconAVGNEG = [];
                        armconAVGNEG = [];
                        
                        estANEG(:,k) = estA(:,k) * (pdfdirection_index(k) > DIRthresh);   
                        estAPOS(:,k) = estA(:,k) * (pdfdirection_index(k) < -DIRthresh);  
                                                
                        RconPOS=estAPOS(Rsites,k);
                        armconPOS=estAPOS(armsites,k);
                        RconAVGPOS = (mean(estAPOS(Rsites,k)));
                        armconAVGPOS = (mean(estAPOS(armsites,k)));
                        RsiteBIASPOS = [RsiteBIASPOS (RconAVGPOS-armconAVGPOS)];
                        
                        RconNEG=estANEG(Rsites,k);
                        armconNEG=estANEG(armsites,k);
                        RconAVGNEG = (mean(estANEG(Rsites,k)));
                        armconAVGNEG = (mean(estANEG(armsites,k)));
                        RsiteBIASNEG = [RsiteBIASNEG (RconAVGNEG-armconAVGNEG)];
                    end
                    
                    %step 8: across replay bins, acquire replay bin spatial content associated with VTA unit spikes
                    
                    for k = 1:length(replaybins) 
                        Vunitpos_nothreshPOS = zeros(size(bins10,2),size(VTAunits,2)); 
                        Vunitpos_nothreshNEG = zeros(size(bins10,2),size(VTAunits,2)); 
                        VunitRconAVGPOS = [];
                        VunitarmconAVGPOS = [];
                        VunitRconAVGNEG = [];
                        VunitarmconAVGNEG = [];
                        
                        x10 = binT10(LEstartime(i))+(k-1)*BINSIZE + delay;
                        y10 = binT10(LEstartime(i))+(k*BINSIZE) + delay;
                        [maxest,Imaxest] = max(estA(:,k));
                        
                        %across VTA units
                        for z=1:size(VTAunits,2) 
                            weight =(length(find(VTAunits(z).spikes>=x10,1,'first'):find(VTAunits(z).spikes<=y10,1,'last'))); 
                            if weight>0
                                Vunitpos_nothreshPOS(:,z) = Vunitpos_nothreshPOS(:,z) + ...
                                    estAPOS(:,k)*weight;
                                
                                Vunitpos_nothreshNEG(:,z) = Vunitpos_nothreshNEG(:,z) + ...
                                    estANEG(:,k)*weight;
                            end
                        end
                        
                        VunitRconAVGPOS = (mean(Vunitpos_nothreshPOS(Rsites,:)));
                        VunitarmconAVGPOS = (mean(Vunitpos_nothreshPOS(armsites,:)));
                        VunitRconAVGPOS(VunitRconAVGPOS==0)=nan;                                                
                        VunitarmconAVGPOS(VunitarmconAVGPOS==0)=nan;
                        VunitRsiteBIASPOS = [VunitRsiteBIASPOS (VunitRconAVGPOS-VunitarmconAVGPOS)'];            
                        
                        VunitRconAVGNEG = (mean(Vunitpos_nothreshNEG(Rsites,:)));
                        VunitarmconAVGNEG = (mean(Vunitpos_nothreshNEG(armsites,:)));
                        VunitRconAVGNEG(VunitRconAVGNEG==0)=nan;
                        VunitarmconAVGNEG(VunitarmconAVGNEG==0)=nan;
                        VunitRsiteBIASNEG = [VunitRsiteBIASNEG (VunitRconAVGNEG-VunitarmconAVGNEG)'];        
                    end
                end
            end
        end
    end
  
elseif strcmp(select_replayslope,'incoming') % select incoming replay events  
    for i=1:length(SPWReventstartime)-1;  % across all SPWR events
        estA = [];
        estD = [];
        
        %step 9: % select replay events occurring at the specified location
        %(Force reward locations), to derive pdf's of position and direction
        
        if intersect(SPWReventstartime(i), replayeventstartime) 
            Vunitpos_nothresh = zeros(size(bins10,2),size(VTAunits,2));
            xi = sum(candidatebinlength(1:i))+1;
            yi = sum(candidatebinlength(1:i+1));
            replaybins = sum(candidatebinlength(1:i))+1:sum(candidatebinlength(1:i+1));
            
            %single replay event spatial pdf x time 
            estA = pdfall10(bins10,xi:yi);  
            
            %single replay event direction pdf x time 
            estD = pdfdirection10(:,xi:yi); 
            
            %inbound direction estimate
            inbound = estD(1,:);  
            
            %outbound direction estimate
            outbound = estD(2,:); 
            
            pdfdirection_index = (outbound-inbound)./(outbound+inbound);
            
            estAPOS = zeros(size(estA));
            estANEG = zeros(size(estA));
            
            if ~strcmp(select_tracktype, 'linear') %consider SWM task
                
                %for incoming replay events with rat at force reward sites 
                %(outbound direction estimate = reverse replay, inbound direction estimate = forward replay)
                
                if replayslope(i) < 0 
                    for k = 1:length(replaybins)
                        
                        %step 10: %across replay bins, compute POS and NEG reward site bias in each bin
                        
                        RconAVGPOS = [];
                        armconAVGPOS = [];
                        
                        RconAVGNEG = [];
                        armconAVGNEG = [];
                        
                        estAPOS(:,k) = estA(:,k) * (pdfdirection_index(k) > DIRthresh);   
                        estANEG(:,k) = estA(:,k) * (pdfdirection_index(k) < -DIRthresh);  
                        
                        RconAVGPOS = (mean(estAPOS(Rsites,k)));
                        armconAVGPOS = (mean(estAPOS(armsites,k)));
                        RsiteBIASPOS = [RsiteBIASPOS (RconAVGPOS-armconAVGPOS)]; 
                        
                        RconAVGNEG = (mean(estANEG(Rsites,k)));
                        armconAVGNEG = (mean(estANEG(armsites,k)));
                        RsiteBIASNEG = [RsiteBIASNEG (RconAVGNEG-armconAVGNEG)]; 
                    end
                    
                    %step 11: across replay bins and VTA units, acquire replay bin spatial content associated with VTA unit spikes
                    
                    for k = 1:length(replaybins) 
                        Vunitpos_nothreshPOS = zeros(size(bins10,2),size(VTAunits,2)); 
                        Vunitpos_nothreshNEG = zeros(size(bins10,2),size(VTAunits,2)); 
                        
                        VunitRconAVGPOS = [];
                        VunitarmconAVGPOS = [];
                        
                        VunitRconAVGNEG = [];
                        VunitarmconAVGNEG = [];
                        
                        x10 = binT10(LEstartime(i))+(k-1)*BINSIZE + delay;
                        y10 = binT10(LEstartime(i))+(k*BINSIZE) + delay;
                        [maxest,Imaxest] = max(estA(:,k));
                        
                        %across VTA units
                        for z=1:size(VTAunits,2)                           
                            weight =(length(find(VTAunits(z).spikes>=x10,1,'first'):find(VTAunits(z).spikes<=y10,1,'last'))); 
                            if weight>0
                                Vunitpos_nothreshPOS(:,z) = Vunitpos_nothreshPOS(:,z) + ...
                                    estAPOS(:,k)*weight;
                                
                                Vunitpos_nothreshNEG(:,z) = Vunitpos_nothreshNEG(:,z) + ...
                                    estANEG(:,k)*weight;
                            end
                        end
                        
                        VunitRconAVGPOS = (mean(Vunitpos_nothreshPOS(Rsites,:)));
                        VunitarmconAVGPOS = (mean(Vunitpos_nothreshPOS(armsites,:)));
                        VunitRconAVGPOS(VunitRconAVGPOS==0)=nan;                                                
                        VunitarmconAVGPOS(VunitarmconAVGPOS==0)=nan;                                            
                        VunitRsiteBIASPOS = [VunitRsiteBIASPOS (VunitRconAVGPOS-VunitarmconAVGPOS)'];           
                        
                        VunitRconAVGNEG = (mean(Vunitpos_nothreshNEG(Rsites,:)));
                        VunitarmconAVGNEG = (mean(Vunitpos_nothreshNEG(armsites,:)));
                        VunitRconAVGNEG(VunitRconAVGNEG==0)=nan;                                                
                        VunitarmconAVGNEG(VunitarmconAVGNEG==0)=nan;                                            
                        VunitRsiteBIASNEG = [VunitRsiteBIASNEG (VunitRconAVGNEG-VunitarmconAVGNEG)'];           
                    end
                end  
                
            elseif strcmp(select_tracktype, 'linear') %consider linear track
                
                %for incoming replay events with rat at near (left) reward site.
                %(inbound = forward replay; outbound = reverse replay)
                
                if replayslope(i) < 0 && positionatreplay(i) < bins(fix(length(bins)/2)) 
                    
                    %step 12: across replay bins, compute POS and NEG reward site bias in each bin
                    
                    for k = 1:length(replaybins) 
                        RconAVGPOS = [];
                        armconAVGPOS = [];
                        
                        RconAVGNEG = [];
                        armconAVGNEG = [];
                        
                        estAPOS(:,k) = estA(:,k) * (pdfdirection_index(k) > DIRthresh);  
                        estANEG(:,k) = estA(:,k) * (pdfdirection_index(k) < -DIRthresh); 
                        
                        RconAVGPOS = (mean(estAPOS(Rsites,k)));
                        armconAVGPOS = (mean(estAPOS(armsites,k)));
                        RsiteBIASPOS = [RsiteBIASPOS (RconAVGPOS-armconAVGPOS)]; 
                        
                        RconAVGNEG = (mean(estANEG(Rsites,k)));
                        armconAVGNEG = (mean(estANEG(armsites,k)));
                        RsiteBIASNEG = [RsiteBIASNEG (RconAVGNEG-armconAVGNEG)]; 
                    end
                    
                    %step 13: across replay bins and VTA units, acquire replay bin spatial content associated with VTA unit spikes
                    
                    for k = 1:length(replaybins) 
                        Vunitpos_nothreshPOS = zeros(size(bins10,2),size(VTAunits,2)); 
                        Vunitpos_nothreshNEG = zeros(size(bins10,2),size(VTAunits,2)); 
                        VunitRconAVGPOS = [];
                        VunitarmconAVGPOS = [];
                        
                        VunitRconAVGNEG = [];
                        VunitarmconAVGNEG = [];
                        
                        x10 = binT10(LEstartime(i))+(k-1)*BINSIZE + delay;
                        y10 = binT10(LEstartime(i))+(k*BINSIZE) + delay;
                        [maxest,Imaxest] = max(estA(:,k));
                        
                        %across VTA units
                        for z=1:size(VTAunits,2) 
                            weight =(length(find(VTAunits(z).spikes>=x10,1,'first'):find(VTAunits(z).spikes<=y10,1,'last')));  
                            if weight>0
                                Vunitpos_nothreshPOS(:,z) = Vunitpos_nothreshPOS(:,z) + ...
                                    estAPOS(:,k)*weight;
                                
                                Vunitpos_nothreshNEG(:,z) = Vunitpos_nothreshNEG(:,z) + ...
                                    estANEG(:,k)*weight;
                            end
                        end
                        
                        VunitRconAVGPOS = (mean(Vunitpos_nothreshPOS(Rsites,:)));
                        VunitarmconAVGPOS = (mean(Vunitpos_nothreshPOS(armsites,:)));
                        VunitRconAVGPOS(VunitRconAVGPOS==0)=nan;
                        VunitarmconAVGPOS(VunitarmconAVGPOS==0)=nan;
                        VunitRsiteBIASPOS = [VunitRsiteBIASPOS (VunitRconAVGPOS-VunitarmconAVGPOS)'];            
                        
                        VunitRconAVGNEG = (mean(Vunitpos_nothreshNEG(Rsites,:)));
                        VunitarmconAVGNEG = (mean(Vunitpos_nothreshNEG(armsites,:)));
                        VunitRconAVGNEG(VunitRconAVGNEG==0)=nan;
                        VunitarmconAVGNEG(VunitarmconAVGNEG==0)=nan;
                        VunitRsiteBIASNEG = [VunitRsiteBIASNEG (VunitRconAVGNEG-VunitarmconAVGNEG)'];            
                    end
                    
                %consider the case of outgoing decoded replay events with rat at far (right) reward site. 
                % (Outbound = forward replay; Inbound = reverse replay)
                
                elseif replayslope(i) > 0 && positionatreplay(i) > bins(fix(length(bins)/2)) 
                    
                    %step 14: across replay bins, compute POS and NEG reward site bias in each bin
                    
                    for k = 1:length(replaybins) 
                        RconAVGPOS = [];
                        armconAVGPOS = [];
                        
                        RconAVGNEG = [];
                        armconAVGNEG = [];
                        
                        estANEG(:,k) = estA(:,k) * (pdfdirection_index(k) > DIRthresh);   
                        estAPOS(:,k) = estA(:,k) * (pdfdirection_index(k) < -DIRthresh);  
                                                
                        RconPOS=estAPOS(Rsites,k);
                        armconPOS=estAPOS(armsites,k);
                        RconAVGPOS = (mean(estAPOS(Rsites,k)));
                        armconAVGPOS = (mean(estAPOS(armsites,k)));
                        RsiteBIASPOS = [RsiteBIASPOS (RconAVGPOS-armconAVGPOS)];
                        
                        RconNEG=estANEG(Rsites,k);
                        armconNEG=estANEG(armsites,k);
                        RconAVGNEG = (mean(estANEG(Rsites,k)));
                        armconAVGNEG = (mean(estANEG(armsites,k)));
                        RsiteBIASNEG = [RsiteBIASNEG (RconAVGNEG-armconAVGNEG)];
                    end
                    
                    %step 15: across replay bins units, acquire replay bin spatial content associated with VTA unit spikes
                    
                    for k = 1:length(replaybins) 
                        Vunitpos_nothreshPOS = zeros(size(bins10,2),size(VTAunits,2)); 
                        Vunitpos_nothreshNEG = zeros(size(bins10,2),size(VTAunits,2)); 
                        VunitRconAVGPOS = [];
                        VunitarmconAVGPOS = [];
                        VunitRconAVGNEG = [];
                        VunitarmconAVGNEG = [];
                        
                        x10 = binT10(LEstartime(i))+(k-1)*BINSIZE + delay;
                        y10 = binT10(LEstartime(i))+(k*BINSIZE) + delay;
                        [maxest,Imaxest] = max(estA(:,k));
                        
                        for z=1:size(VTAunits,2) %across VTA units
                            weight =(length(find(VTAunits(z).spikes>=x10,1,'first'):find(VTAunits(z).spikes<=y10,1,'last'))); 
                            if weight>0
                                Vunitpos_nothreshPOS(:,z) = Vunitpos_nothreshPOS(:,z) + ...
                                    estAPOS(:,k)*weight;
                                
                                Vunitpos_nothreshNEG(:,z) = Vunitpos_nothreshNEG(:,z) + ...
                                    estANEG(:,k)*weight;
                            end
                        end
                        
                        VunitRconAVGPOS = (mean(Vunitpos_nothreshPOS(Rsites,:)));
                        VunitarmconAVGPOS = (mean(Vunitpos_nothreshPOS(armsites,:)));
                        VunitRconAVGPOS(VunitRconAVGPOS==0)=nan;                                                
                        VunitarmconAVGPOS(VunitarmconAVGPOS==0)=nan;
                        VunitRsiteBIASPOS = [VunitRsiteBIASPOS (VunitRconAVGPOS-VunitarmconAVGPOS)'];            
                        
                        VunitRconAVGNEG = (mean(Vunitpos_nothreshNEG(Rsites,:)));
                        VunitarmconAVGNEG = (mean(Vunitpos_nothreshNEG(armsites,:)));
                        VunitRconAVGNEG(VunitRconAVGNEG==0)=nan;
                        VunitarmconAVGNEG(VunitarmconAVGNEG==0)=nan;
                        VunitRsiteBIASNEG = [VunitRsiteBIASNEG (VunitRconAVGNEG-VunitarmconAVGNEG)'];            
                    end
                end
            end
        end
    end
end

%step 16: transform Rsitebias and VunitRsitebias into binary output

RsiteBIASPOS(RsiteBIASPOS==0) = nan;
RsiteBIASPOS(RsiteBIASPOS>0) = 1;
RsiteBIASPOS(RsiteBIASPOS<0) = 0;

RsiteBIASNEG(RsiteBIASNEG==0) = nan;
RsiteBIASNEG(RsiteBIASNEG>0) = 1;
RsiteBIASNEG(RsiteBIASNEG<0) = 0;

VunitRsiteBIASPOS(VunitRsiteBIASPOS>0) = 1;
VunitRsiteBIASPOS(VunitRsiteBIASPOS<0) = 0;

VunitRsiteBIASNEG(VunitRsiteBIASNEG>0) = 1;
VunitRsiteBIASNEG(VunitRsiteBIASNEG<0) = 0;

%step 17: assign RsiteBIAS and VunitRsiteBIAS

if strcmp(select_direction, 'POS')
    RsiteBIAS=RsiteBIASPOS;
    VunitRsiteBIAS=VunitRsiteBIASPOS;
elseif strcmp(select_direction, 'NEG')
    RsiteBIAS=RsiteBIASNEG;
    VunitRsiteBIAS=VunitRsiteBIASNEG;
end

%step 18: after computing RsiteBIAS and VunitRsiteBIAS for each of the 4 permutations 
% (1. POS, outgoing; 2. POS, incoming; 3. NEG, outgoing; 4. NEG, incoming), 
% derive RsiteBIAS and VunitRsiteBIAS for forward/reverse centrifugal/centripetal replay as follows:

% forward_replay  { select_direction = 'POS', select_replayslope = 'outgoing'; 
%                   select_direction = 'NEG', select_replayslope = 'incoming'}

RsiteBIASPOSfwd     = [Rsitebias1 Rsitebias4];
VunitRsiteBIASfwd   = [VunitRsiteBIAS1 VunitRsiteBIAS4];

% reverse_replay  { select_direction = 'NEG', select_replayslope = 'outgoing'; 
%                   select_direction = 'POS', select_replayslope = 'incoming'}

RsiteBIASPOSrvs     = [Rsitebias2 Rsitebias3];
VunitRsiteBIASrvs   = [VunitRsiteBIAS2 VunitRsiteBIAS3];


% centrifugal_replay  { select_direction = 'POS', select_replayslope = 'outgoing'; 
%                       select_direction = 'NEG', select_replayslope = 'outgoing'}

RsiteBIASPOScentrifugal     = [Rsitebias1 Rsitebias3];
VunitRsiteBIAScentrifugal   = [VunitRsiteBIAS1 VunitRsiteBIAS3];

% centripetal_replay  { select_direction = 'NEG', select_replayslope = 'incoming'; 
%                       select_direction = 'POS', select_replayslope = 'incoming'}

RsiteBIASPOScentripetal     = [Rsitebias2 Rsitebias4];
VunitRsiteBIAScentripetal   = [VunitRsiteBIAS2 VunitRsiteBIAS4];


