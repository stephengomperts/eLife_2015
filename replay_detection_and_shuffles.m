% replay detection and shuffles to find significant replay events

% This code computes the probability that linear spatial trajectories of 
% replayed position during SPW-R events are significantly different from chance. 
% Four trajectories are evaluated, given the four trajectories possible on 
% the SWM task, as defined in Supplementary Figure 3a: ABCD, ABCF, EBCD, EBCF.
% It outputs the indices of SPW-R events that are replay events and that are nonreplay events.
% It uses a line finding algorithm (implemented below as est_line_detect) based on Davidson et al. 2009 and further 
% described in Kloosterman, 2012 that uses the radon transform to find the best 
% line in the spatial pdf of each candidate replay event.

%prerequisite data:
% SPWReventstartime :      vector of SPW-R event start times (sec)
% startime :               start time of session (sec)
% stoptime :               stop time of session (sec)
% BINSIZE :                replay event time bin size (25 msec)
% pdfall :                 array of spatially decoded, concatenated candidate replay events, comprising a spatial pdf (10cm bins) x time (25msec bins) 
% candidatebinlength :     vector of SPW-R time bin lengths

% bins10 :                 vector of the 10cm spatial bins that span the track
% binsABCD :               vector of spatial bins for trajectory ABCD
% edgeABCD :               vector of spatial nodes for trajectory ABCD
% binsABCF :               vector of spatial bins for trajectory ABCF
% edgeABCF :               vector of spatial nodes for trajectory ABCF
% binsEBCD :               vector of spatial bins for trajectory EBCD
% edgeEBCD :               vector of spatial nodes for trajectory EBCD
% binsEBCF :               vector of spatial bins for trajectory EBCF
% edgeEBCF :               vector of spatial nodes for trajectory EBCF


maxcolumn = 1500;                       %define the number of column shuffles
maxpseudo = 1500;                       %define the number of pseudoevent shuffles
numevents = length(SPWReventstartime);  

score = [];

slopeALLABCD     = [];
interceptALLABCD = [];
scoreALLABCD     = [];
slopeA     = [];
interceptA = [];
scoreA     = [];

slopeALLABCDrand     = [];
interceptALLABCDrand = [];
scoreALLABCDrand     = [];
slopeAr     = [];
interceptAr = [];
scoreAr     = [];

scoreALLABCDpseudorand = [];
interceptALLABCDpseudorand = [];
slopeApr     = [];
interceptApr = [];
scoreApr     = [];

slopeALLABCF     = [];
interceptALLABCF = [];
scoreALLABCF     = [];
slopeALLABCFrand     = [];
interceptALLABCFrand = [];
scoreALLABCFrand     = [];
interceptALLABCFpseudorand = [];
scoreALLABCFpseudorand     = [];

slopeALLEBCD     = [];
interceptALLEBCD = [];
scoreALLEBCD     = [];
slopeALLEBCDrand     = [];
interceptALLEBCDrand = [];
scoreALLEBCDrand     = [];
interceptALLEBCDpseudorand = [];
scoreALLEBCDpseudorand     = [];

slopeALLEBCF     = [];
interceptALLEBCF = [];
scoreALLEBCF     = [];
slopeALLEBCFrand     = [];
interceptALLEBCFrand = [];
scoreALLEBCFrand     = [];
interceptALLEBCFpseudorand = [];
scoreALLEBCFpseudorand     = [];

appendedestA = pdfall10;
appendedestAABCD = appendedestA(binsABCD,:);
appendedestAABCF = appendedestA(binsABCF,:);
appendedestAEBCD = appendedestA(binsEBCD,:);
appendedestAEBCF = appendedestA(binsEBCF,:);

scoreindex = [];
slopebest     = [];
interceptbest = [];
scorebest    = [];
maxscore = [];
maxscoreindex = [];
for i = 1:numevents %loop through SPW-R events
    
    %step 1: compute replay line scores for each trajectory
    
    disp(i)
    xi = sum(candidatebinlength(1:i))+1;
    yi = sum(candidatebinlength(1:i+1));
    replaybins = sum(candidatebinlength(1:i))+1:sum(candidatebinlength(1:i+1)); % indices of pdfall for candidate replay event i
    
    pdfstartindex = xi;
    pdfendindex   = yi;
    
    estAABCD = pdfall10(binsABCD,pdfstartindex:pdfendindex);
    [slopeA,interceptA,scoreA] = est_line_detect(replaybins,binsABCD,estAABCD, 'kernelwidth',[3 0]);
    slopeALLABCD = [slopeALLABCD; slopeA];
    interceptALLABCD = [interceptALLABCD; interceptA];
    scoreALLABCD = [scoreALLABCD; scoreA]; %replay line score
    
    estAABCF = pdfall10(binsABCF,pdfstartindex:pdfendindex);
    [slopeA,interceptA,scoreA] = est_line_detect(replaybins,binsABCF,estAABCF, 'kernelwidth',[3 0]);
    slopeALLABCF = [slopeALLABCF; slopeA];
    interceptALLABCF = [interceptALLABCF; interceptA];
    scoreALLABCF = [scoreALLABCF; scoreA]; %replay line score
    
    estAEBCD = pdfall10(binsEBCD,pdfstartindex:pdfendindex);
    [slopeA,interceptA,scoreA] = est_line_detect(replaybins,binsEBCD,estAEBCD, 'kernelwidth',[3 0]);
    slopeALLEBCD = [slopeALLEBCD; slopeA];
    interceptALLEBCD = [interceptALLEBCD; interceptA];
    scoreALLEBCD = [scoreALLEBCD; scoreA]; %replay line score
    
    estAEBCF = pdfall10(binsEBCF,pdfstartindex:pdfendindex);
    [slopeA,interceptA,scoreA] = est_line_detect(replaybins,binsEBCF,estAEBCF, 'kernelwidth',[3 0]);
    slopeALLEBCF = [slopeALLEBCF; slopeA];
    interceptALLEBCF = [interceptALLEBCF; interceptA];
    scoreALLEBCF = [scoreALLEBCF; scoreA]; %replay line score
    
    %step 2: keep the best replay line score
    
    [maxscore,maxscoreindex] = max([scoreALLABCD(i) scoreALLABCF(i) scoreALLEBCD(i),scoreALLEBCF(i)]);  
    scoreindex = [scoreindex; maxscoreindex];
    if maxscoreindex == 1
        slopebest = [slopebest; slopeALLABCD(i)];
        interceptbest = [interceptbest; interceptALLABCD(i)];
        scorebest = [scorebest; scoreALLABCD(i)];
    elseif maxscoreindex == 2
        slopebest = [slopebest; slopeALLABCF(i)];
        interceptbest = [interceptbest; interceptALLABCF(i)];
        scorebest = [scorebest; scoreALLABCF(i)];
    elseif maxscoreindex == 3
        slopebest = [slopebest; slopeALLEBCD(i)];
        interceptbest = [interceptbest; interceptALLEBCD(i)];
        scorebest = [scorebest; scoreALLEBCD(i)];
    elseif maxscoreindex == 4
        slopebest = [slopebest; slopeALLEBCF(i)];
        interceptbest = [interceptbest; interceptALLEBCF(i)];
        scorebest = [scorebest; scoreALLEBCF(i)];
    end
    
    %step 3. column shuffle
    
    scorerandindex(i).event = [];
    slopebestrand(i).event = [];
    interceptbestrand(i).event = [];
    scorebestrand(i).event = [];
    for j=1:maxcolumn
        estAABCDrand = randomize(estAABCD,'method','cycle');
        [slopeAr,interceptAr,scoreAr] = est_line_detect(replaybins,binsABCD,estAABCDrand, 'kernelwidth',[3 0]);
        slopeALLABCDrand(i).event(j)     = slopeAr;
        interceptALLABCDrand(i).event(j) = interceptAr;
        scoreALLABCDrand(i).event(j)     = scoreAr;
        
        estAABCFrand = randomize(estAABCF,'method','cycle');
        [slopeAr,interceptAr,scoreAr] = est_line_detect(replaybins,binsABCF,estAABCFrand, 'kernelwidth',[3 0]);
        slopeALLABCFrand(i).event(j)     = slopeAr;
        interceptALLABCFrand(i).event(j) = interceptAr;
        scoreALLABCFrand(i).event(j)     = scoreAr;
        
        estAEBCDrand = randomize(estAEBCD,'method','cycle');
        [slopeAr,interceptAr,scoreAr] = est_line_detect(replaybins,binsEBCD,estAEBCDrand, 'kernelwidth',[3 0]);
        slopeALLEBCDrand(i).event(j)     = slopeAr;
        interceptALLEBCDrand(i).event(j) = interceptAr;
        scoreALLEBCDrand(i).event(j)     = scoreAr;
        
        estAEBCFrand = randomize(estAEBCF,'method','cycle');
        [slopeAr,interceptAr,scoreAr] = est_line_detect(replaybins,binsEBCF,estAEBCFrand, 'kernelwidth',[3 0]);
        slopeALLEBCFrand(i).event(j)     = slopeAr;
        interceptALLEBCFrand(i).event(j) = interceptAr;
        scoreALLEBCFrand(i).event(j)     = scoreAr;
        
        [maxscorerand,maxscorerandindex] = max([scoreALLABCDrand(i).event(j)   scoreALLABCFrand(i).event(j) ...
            scoreALLEBCDrand(i).event(j)   scoreALLEBCFrand(i).event(j)   ]);
        scorerandindex(i).event = [scorerandindex(i).event; maxscorerandindex];
        if maxscorerandindex == 1
            slopebestrand(i).event = [slopebestrand(i).event; slopeALLABCDrand(i).event(j)];
            interceptbestrand(i).event = [interceptbestrand(i).event; interceptALLABCDrand(i).event(j)];
            scorebestrand(i).event = [scorebestrand(i).event; scoreALLABCDrand(i).event(j)];
        elseif maxscorerandindex == 2
            slopebestrand(i).event = [slopebestrand(i).event; slopeALLABCFrand(i).event(j)];
            interceptbestrand(i).event = [interceptbestrand(i).event; interceptALLABCFrand(i).event(j)];
            scorebestrand(i).event = [scorebestrand(i).event; scoreALLABCFrand(i).event(j)];
        elseif maxscorerandindex == 3
            slopebestrand(i).event = [slopebestrand(i).event; slopeALLEBCDrand(i).event(j)];
            interceptbestrand(i).event = [interceptbestrand(i).event; interceptALLEBCDrand(i).event(j)];
            scorebestrand(i).event = [scorebestrand(i).event; scoreALLEBCDrand(i).event(j)];
        elseif maxscorerandindex == 4
            slopebestrand(i).event = [slopebestrand(i).event; slopeALLEBCFrand(i).event(j)];
            interceptbestrand(i).event = [interceptbestrand(i).event; interceptALLEBCFrand(i).event(j)];
            scorebestrand(i).event = [scorebestrand(i).event; scoreALLEBCFrand(i).event(j)];
        end
               
    end
    
    %step 4. compute probability of best linescore relative to column shuffle
    
    scorebins = linspace(0,max(scorebestrand(i).event)+0.02,100);
    dist=histc(scorebestrand(i).event,scorebins);
    score.columnshuffle(i,1).pBEST = numel(find(scorebestrand(i).event>scorebest(i)))/maxcolumn; % probability of best linescore (normalized by number of column shuffles)
    score.columnshuffle(i).distBEST = scorebestrand(i).event';
    score.columnshuffle(i).bindistBEST = dist';
    score.columnshuffle(i).binsBEST = scorebins';
  
    %step 5. pseudoevent shuffle
    
    scorepseudorandindex(i).event = [];
    slopebestpseudorand(i).event = [];
    interceptbestpseudorand(i).event = [];
    scorebestpseudorand(i).event = [];
    for j=1:maxpseudo
        k = length(replaybins);
        n = length(appendedestAABCD(1,:));
        y1 = randsample(n,k);
        estAABCDpseudo = appendedestAABCD(:,y1);
        [slopeALLABCDpseudorand(i).event(j),interceptALLABCDpseudorand(i).event(j),scoreApr] = est_line_detect(replaybins,binsABCD,estAABCDpseudo, 'kernelwidth',[3 0]);
        scoreALLABCDpseudorand(i).event(j)     = scoreApr;
        
        n = length(appendedestAABCF(1,:));
        y1 = randsample(n,k);
        estAABCFpseudo = appendedestAABCF(:,y1);
        [slopeALLABCFpseudorand(i).event(j),interceptALLABCFpseudorand(i).event(j),scoreApr] = est_line_detect(replaybins,binsABCF,estAABCFpseudo, 'kernelwidth',[3 0]);
        scoreALLABCFpseudorand(i).event(j)     = scoreApr;
        
        n = length(appendedestAEBCD(1,:));
        y1 = randsample(n,k);
        estAEBCDpseudo = appendedestAEBCD(:,y1);
        [slopeALLEBCDpseudorand(i).event(j),interceptALLEBCDpseudorand(i).event(j),scoreApr] = est_line_detect(replaybins,binsEBCD,estAEBCDpseudo, 'kernelwidth',[3 0]);
        scoreALLEBCDpseudorand(i).event(j)     = scoreApr;
        
        n = length(appendedestAEBCF(1,:));
        y1 = randsample(n,k);
        estAEBCFpseudo = appendedestAEBCF(:,y1);
        [slopeALLEBCFpseudorand(i).event(j),interceptALLEBCFpseudorand(i).event(j),scoreApr] = est_line_detect(replaybins,binsEBCF,estAEBCFpseudo, 'kernelwidth',[3 0]);
        scoreALLEBCFpseudorand(i).event(j)     = scoreApr;
        
         [maxscorepseudorand,maxscorepseudorandindex] = max([scoreALLABCDpseudorand(i).event(j)   scoreALLABCFpseudorand(i).event(j) ...
            scoreALLEBCDpseudorand(i).event(j)   scoreALLEBCFpseudorand(i).event(j)   ]);
        scorepseudorandindex(i).event = [scorepseudorandindex(i).event; maxscorepseudorandindex];
        if maxscorepseudorandindex == 1
            slopebestpseudorand(i).event = [slopebestpseudorand(i).event; slopeALLABCDpseudorand(i).event(j)];
            interceptbestpseudorand(i).event = [interceptbestpseudorand(i).event; interceptALLABCDpseudorand(i).event(j)];
            scorebestpseudorand(i).event = [scorebestpseudorand(i).event; scoreALLABCDpseudorand(i).event(j)];
        elseif maxscorepseudorandindex == 2
            slopebestpseudorand(i).event = [slopebestpseudorand(i).event; slopeALLABCFpseudorand(i).event(j)];
            interceptbestpseudorand(i).event = [interceptbestpseudorand(i).event; interceptALLABCFpseudorand(i).event(j)];
            scorebestpseudorand(i).event = [scorebestpseudorand(i).event; scoreALLABCFpseudorand(i).event(j)];
        elseif maxscorepseudorandindex == 3
            slopebestpseudorand(i).event = [slopebestpseudorand(i).event; slopeALLEBCDpseudorand(i).event(j)];
            interceptbestpseudorand(i).event = [interceptbestpseudorand(i).event; interceptALLEBCDpseudorand(i).event(j)];
            scorebestpseudorand(i).event = [scorebestpseudorand(i).event; scoreALLEBCDpseudorand(i).event(j)];
        elseif maxscorepseudorandindex == 4
            slopebestpseudorand(i).event = [slopebestpseudorand(i).event; slopeALLEBCFpseudorand(i).event(j)];
            interceptbestpseudorand(i).event = [interceptbestpseudorand(i).event; interceptALLEBCFpseudorand(i).event(j)];
            scorebestpseudorand(i).event = [scorebestpseudorand(i).event; scoreALLEBCFpseudorand(i).event(j)];
        end       
        
        
    end
    
    %step 6. compute probability of best linescore relative to pseudoevent  shuffle
    
    scorebins = linspace(0,max(scorebestpseudorand(i).event)+0.02,100);
    dist=histc(scorebestpseudorand(i).event,scorebins);
    score.pseudoshuffle(i,1).pBEST = numel(find(scorebestpseudorand(i).event>scorebest(i)))/maxpseudo; % probability of best linescore (normalized by number of pseudoevent shuffles)
    score.pseudoshuffle(i).distBEST = scorebestpseudorand(i).event';
    score.pseudoshuffle(i).bindistBEST = dist';
    score.pseudoshuffle(i).binsBEST = scorebins';
    
    score.pBEST.columnshuffle(i) = score.columnshuffle(i,1).pBEST;    
    score.pBEST.pseudoshuffle(i) = score.pseudoshuffle(i,1).pBEST;    
        
end

%step 7. find significant replay trajectories

pthresh = 0.05;
coBEST=find(score.pBEST.columnshuffle < pthresh);
psBEST=find(score.pBEST.pseudoshuffle < pthresh);
replayevent_idx = intersect(coBEST,psBEST);         %indices of SPW-R events associated with significant (at p<0.05 level) replay of spatial trajectories

%step 8. find nonsignificant replay trajectories

coNON=find(score.pBEST.columnshuffle > 0.2);
psNON=find(score.pBEST.pseudoshuffle > 0.2);
nonreplayevent_idx = intersect(coNON,psNON);        %indices of SPW-R events associated with nonsignificant replay of spatial trajectories



