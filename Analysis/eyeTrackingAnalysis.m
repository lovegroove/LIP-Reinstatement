% Analyze eye tracking data
% PDS.eyepos = [x,y,trialTime, trialState], each cell from a trial
% PDS.data.eyeposloop is the same i think
% dv.states.SHOWPROBE = 11;
% objLocs = % [correctObj, foilObj, loc1, loc2, correctLoc]
% PDS.data.correctObject(iTrial)
   
    format longg
    
    eyePos = PDS.eyepos(strcmp(PDS.trialType,'test'));
    objLocs = PDS.data.objectLocs(strcmp(PDS.trialType,'test'));
    
    state = dv.states.SHOWPROBE; %only the eyePos during memory probe (just in case we want to switch what state we are analyzing, BUT the areas where the shapes are will be different in other states of the trial)
    
    % NOTE: zeros in cells in eyePosProp are from break fixation, there are no eye positions during the show probe state then
    
    n = length(eyePos);
    eyePosState = cell(1,n);
    eyePosProp = cell(1,n);
    emptyCells = zeros(n,1);
    firstSaccadeSample = zeros(n,1);
    firstSaccadeTime = zeros(n,1);
    firstSaccadeStatus = zeros(n,1);
    relativeMatchTime = zeros(n,1);
    propMatchTime = zeros(n,1);
    propFoilTime = zeros(n,1);
    propNeitherTime = zeros(n,1);
    relativeRatio = zeros(n,1);
    
    for iTrial = 1:n
        
        eyePosState{iTrial} = eyePos{iTrial}(eyePos{iTrial}(:,4)==state,:); % only eye position during a certain state
        sampleLen = length(eyePosState{iTrial});
        
        for jSample = 1:sampleLen
            % eyes within 1st object boundaries
            if eyePosState{iTrial}(jSample,1) >= objLocs{iTrial}{1,3}(1,1) && eyePosState{iTrial}(jSample,1) <= objLocs{iTrial}{1,3}(1,3) && eyePosState{iTrial}(jSample,2) >= objLocs{iTrial}{1,3}(1,2) && eyePosState{iTrial}(jSample,2) <= objLocs{iTrial}{1,3}(1,4)
                
                if objLocs{iTrial}{1,5} == 1
                    % correct
                    eyePosProp{iTrial}(jSample,1) = 1;
                    
                elseif objLocs{iTrial}{1,5} == 2
                    % incorrect (FA)
                    eyePosProp{iTrial}(jSample,1) = 2;
                end
                
                % eyes within 2nd object boundaries
            elseif eyePosState{iTrial}(jSample,1) >= objLocs{iTrial}{1,4}(1,1) && eyePosState{iTrial}(jSample,1) <= objLocs{iTrial}{1,4}(1,3) && eyePosState{iTrial}(jSample,2) >= objLocs{iTrial}{1,4}(1,2) && eyePosState{iTrial}(jSample,2) <= objLocs{iTrial}{1,4}(1,4)
                
                if objLocs{iTrial}{1,5} == 2
                    % correct
                    eyePosProp{iTrial}(jSample,1) = 1;
                    
                elseif objLocs{iTrial}{1,5} == 1
                    % incorrect (FA)
                    eyePosProp{iTrial}(jSample,1) = 2;
                end
                
            else  % eyes in neither
                eyePosProp{iTrial}(jSample,1) = 0;
            end
        end
        
        % flag empty cells due to break fixation
        emptyCells(iTrial) = isempty(eyePosProp{iTrial}); % don't really need anymore, is this still useful?
        
        
        if all(eyePosProp{iTrial}(:)==0) % to deal with weird all zero cell
            eyePosProp{iTrial} = [];
        end
        
        if ~isempty(eyePosProp{iTrial}) 
            
            % first saccade
            firstSaccadeSample(iTrial) = find(eyePosProp{iTrial} ~= 0,1);
            % time of first saccade
            eyePosStateTime = eyePosState{iTrial}(:,3); % time from beginning of trial
            firstSaccadeTime(iTrial) = eyePosStateTime(firstSaccadeSample(iTrial)) - eyePosStateTime(1);
            % status of first saccade (location - target or foil)
            firstSaccadeStatus(iTrial) = eyePosProp{iTrial}(firstSaccadeSample(iTrial));
            
            %%%%%%%%%%%%%%%%%%%%%%
            
            % relative viewing
            sampleTime = (eyePosStateTime(end) - eyePosStateTime(1)) / length(eyePosStateTime);
            propMatchTime(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 1)); % why are these all the same??
            propFoilTime(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 2));
            propNeitherTime(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 0));
            relativeMatchTime(iTrial) = propMatchTime(iTrial) / (propMatchTime(iTrial) + propFoilTime(iTrial));
            relativeRatio(iTrial) = propMatchTime(iTrial) / propFoilTime(iTrial);
        else
            firstSaccadeSample(iTrial) = NaN;
            firstSaccadeTime(iTrial) = NaN;
            firstSaccadeStatus(iTrial) = NaN;
            propMatchTime(iTrial) = NaN;
            propFoilTime(iTrial) = NaN;
            propNeitherTime(iTrial) = NaN;
            relativeMatchTime(iTrial) = NaN;
            relativeRatio(iTrial) = NaN;
        end
        
        
    end
    
    % eliminate NaNs
    firstSaccadeSample = firstSaccadeSample(~isnan(firstSaccadeSample));
    firstSaccadeTime = firstSaccadeTime(~isnan(firstSaccadeTime));
    firstSaccadeStatus = firstSaccadeStatus(~isnan(firstSaccadeStatus));
    propMatchTime = propMatchTime(~isnan(propMatchTime));
    propFoilTime = propFoilTime(~isnan(propFoilTime));
    propNeitherTime = propNeitherTime(~isnan(propNeitherTime));
    relativeMatchTime = relativeMatchTime(~isnan(relativeMatchTime));
    relativeRatio = relativeRatio(~isnan(relativeRatio));

%% Plotting - NOTE I'm currently rewriting ploting variables like x and y
% percent correct
pcTrials = length(find(relativeMatchTime > .5)) / length(relativeMatchTime);

figure;
% Viewing Time

% Total viewing time in a session
y = [propMatchTime propFoilTime propNeitherTime];
subplot(2,2,1), h = barwitherr(std(y),mean(y));
title(sprintf('Total Viewing time - %s', dv.subj))
set(gca,'XTickLabel',{'Match','Foil','Neither'})
set(h,'FaceColor','b');
ylabel('Viewing Time (s)')
subplot(2,2,2), boxplot(y)
title(sprintf('Total Viewing time - %s', dv.subj))
ylabel('Viewing Time (s)')
xtix = {'Match','Foil','Neither'};
xtixloc = [1 2 3];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

%  Relative Viewing TIME match vs. foil
subplot(2,2,3), boxplot(relativeMatchTime)
title(sprintf('Relative Viewing Time: match vs. foil, Percent Correct: %0.5g',pcTrials * 100))
ylabel('Relative Viewing Time: match vs. foil')
subplot(2,2,4), hist(relativeMatchTime)
g = findobj(gca,'Type','patch');
set(g,'FaceColor','b');
title(sprintf('Relative Viewing Time: match vs. foil, Percent Correct: %0.5g',pcTrials * 100))
xlabel('Relative Viewing Time: match vs. foil')
ylabel('Number of Trials')


%% Ratio Index - a lot of divide by zero!!!!
figure;
hist(log(relativeRatio))
title('Ratio of Viewing Time: (=1 match=foil, >1 >match)')
xlabel('(Log) Ratio of Viewing Time: match vs. foil')
ylabel('Number of Trials')

figure;
hist((1 - relativeRatio))
title('Index of Viewing Time')
xlabel('Index of Viewing Time: match vs. foil')
ylabel('Number of Trials')

%% Reaction Times
% RTs Anova
[P,ANOVATAB,STATS] = anova1(firstSaccadeTime, firstSaccadeStatus);
title(sprintf('Reaction Time of First Saccade - Subject: %s, F = %0.5g', dv.subj,ANOVATAB{2,5}))
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ylabel('Reaction Time (s)')


%% logistic regression on RTs and first saccades to Match or Foil 
firstSaccadeStatus(firstSaccadeStatus == 2) = 0; % need a binary variable (after you do this, you can't go back and execute the previous cell)

figure;
[b,dev,stats] = glmfit(firstSaccadeTime,firstSaccadeStatus,'binomial','link','logit');
%yfit = glmval(b,x,'logit');

xx = linspace(0,1.8); % or max(firstSaccadeTime)+max(firstSaccadeTime)*.1
yfit = glmval(b,xx,'logit');
plot(firstSaccadeTime,firstSaccadeStatus,'o',xx,yfit,'-')
title('First Saccades to Match or Foil as a function of Reaction Time')
ylabel('Probability that the First Saccade will be to the Match rather than the Foil')
xlabel('Reaction Time')


%% Per trial time course of relative viewing 

% But what would you look at through the trial - the relative viewing is
% only during the probe and constitutes the 'choice'


%% Pupil Size (or edf.gaze.pupil using the mgl toolbox)

% we are getting this data from the eyelink sample buffer rather than eyepos
% percent signal change (change from average)
% samples
% PDS.data.eyelinkSampleBuffer{1}(12,:)
% -32768 = missing data
% what about zeros


% Pupil Size during test
elSamples = PDS.data.eyelinkSampleBuffer(strcmp(PDS.trialType,'test'));
n = length(elSamples);
pupilSize = cell(1,n);

% which eye - got row from el manual
if strcmp(dv.el.EYE_USED, 'LEFT')
    eyeUsed = 12;
elseif strcmp(dv.el.EYE_USED, 'RIGHT')
    eyeUsed = 13;
else
    error('Neither right eye nor left?')
end

for i = 1:n
    pupilSize{i} = elSamples{i}(eyeUsed,:);
    % Get rid of missing data and zeros
    pupilSize{i} = pupilSize{i}(pupilSize{i} > 0);
    
end

n = length(pupilSize); % new pupilSize length
pupilSizeDev = cell(1,n);
pupilSizeAbsDev = cell(1,n);
pcSigChan = cell(1,n);
pupilSizeZ = cell(1,n);
for i = 1:n
    % deviations from the mean
    pupilSizeDev{i} = pupilSize{i} - mean(pupilSize{i});
    pupilSizeAbsDev{i} = abs(pupilSize{i} - mean(pupilSize{i}));
    
    %pcSigChan{i} = (pupilSize{i} / mean(pupilSize{i})); % how to do this right? * 100 doesn't seem right
    
    % compute Z-score for pupil size (*** makes it run a lot longer)
    psZtrial = zeros(1,length(pupilSize{i}));
    for j = 1:length(pupilSize{i})
        psZtrial(j) = (pupilSize{i}(j) - mean(pupilSize{i})) / std(pupilSize{i});
    end
    pupilSizeZ{i} = psZtrial;
    
end

%% Pupil size during PROBE & CUE
% getting pupil size during probe state
probeStartTimes = PDS.timing.timeShowProbeStart(strcmp(PDS.trialType,'test')); % zeros in cells are break fixation trials that never went to the probe state (remove here and in pupil size data)
probeStartms = cellfun(@(x) x * 1000, probeStartTimes); % convert to ms

% getting pupil size during cue state
cueStartTimes = PDS.timing.timeShowCueStart(strcmp(PDS.trialType,'test')); % zeros in cells are break fixation trials that never went to the probe state (remove here and in pupil size data)
cueStartms = cellfun(@(x) x * 1000, cueStartTimes); % convert to ms

% Getting pupil size during delay state
delayStartTimes = PDS.timing.timeDelayStart(strcmp(PDS.trialType,'test'));
delayStartms = cellfun(@(x) x * 1000, delayStartTimes);

 % remove trials that the subject never made it to the probe state
empInd = cellfun(@isempty, eyePosProp);
chooseIndex = 2; % which index do you want to use to remove trials? 1 = zeros in probe state, 2 = empty cells in eye positions 
switch chooseIndex
    case 1
        pupilSizeDev(:,probeStartms==0) = [];
        pupilSizeAbsDev(:,probeStartms==0) = []; 
        pupilSizeZ(:,probeStartms==0) = []; 
        probeStartms = probeStartms(probeStartms~=0); 
        
        pupilSizeDev(:,cueStartms==0) = [];
        pupilSizeAbsDev(:,cueStartms==0) = []; 
        pupilSizeZ(:,cueStartms==0) = []; 
        cueStartms = cueStartms(cueStartms~=0); 
        
    case 2
        pupilSizeDev(:,empInd==1) = [];
        pupilSizeAbsDev(:,empInd==1) = []; 
        pupilSizeZ(:,empInd==1) = []; 
        probeStartms(empInd==1) = []; 
        cueStartms(empInd==1) = []; 
end

n = length(probeStartms);
psDevProbe = cell(1,n);
psAbsDevProbe = cell(1,n);
psZProbe = cell(1,n);
for i = 1:n
    probeStart = round(probeStartms(i));
    probeEnd = probeStart + (dv.pa.probeTime * 1000); % convert to ms
    psDevProbe{i} = pupilSizeDev{i}(probeStart:probeEnd);
    psAbsDevProbe{i} = pupilSizeAbsDev{i}(probeStart:probeEnd);
    psZProbe{i} = pupilSizeZ{i}(probeStart:probeEnd);
end

n = length(cueStartms);
psDevCue = cell(1,n);
psAbsDevCue = cell(1,n);
psZcue = cell(1,n);
psZcue2probeEnd = cell(1,n);
for i = 1:n
    cueStart = round(cueStartms(i));
    cueEnd = cueStart + (dv.pa.sceneTime * 1000); % convert to ms
    psDevCue{i} = pupilSizeDev{i}(cueStart:cueEnd);
    psAbsDevCue{i} = pupilSizeAbsDev{i}(cueStart:cueEnd);
    psZcue{i} = pupilSizeZ{i}(cueStart:cueEnd);
    psZcue2probeEnd{i} = pupilSizeZ{i}(cueStart:probeEnd); % this might not really get the segment i want
end

% Means
meanPSDevProbe = cellfun(@mean, psDevProbe);
meanPSAbsDevProbe = cellfun(@mean, psAbsDevProbe);
meanPSZprobe = cellfun(@mean, psZProbe);

meanPSDevCue = cellfun(@mean, psDevCue);
meanPSAbsDevCue = cellfun(@mean, psAbsDevCue);
meanPSZcue = cellfun(@mean, psZcue);

%% Pupil Size over the Time Course of the Whole trial
maxSize = max(cellfun(@numel,pupilSizeZ));    %# Get the maximum vector size
aFunc = @(x) [x nan(1,maxSize-numel(x))];  %# Create an anonymous function
psZmat = cellfun(aFunc,pupilSizeZ,'UniformOutput',false);  %# Pad each cell with NaNs
psZmat = vertcat(psZmat{:});                  %# Vertically concatenate cells

meanPSZ = nanmean(psZmat);
stdPSZ = nanstd(psZmat);

figure, shadedErrorBar([],meanPSZ,stdPSZ,'-g',1)
title('Mean Pupil Size')
xlabel('Time (ms)')
ylabel('Pupil Size (Normalized)')
vline(mean(probeStartms), 'b', 'Probe')
vline(mean(cueStartms), 'b', 'Cue')
vline(mean(delayStartms), 'b', 'Delay')

psMeanChoseMatch = nanmean(psZmat(relativeMatchTime > .5,:));
psStdChoseMatch = nanstd(psZmat(relativeMatchTime > .5,:));
psMeanChoseFoil = nanmean(psZmat(relativeMatchTime < .5,:));
psStdChoseFoil = nanstd(psZmat(relativeMatchTime < .5,:));

% Fancey figures - comparing match vs foil
fanceyFigMeans = [psMeanChoseMatch; psMeanChoseFoil];
fanceyFigStds = [psStdChoseMatch; psStdChoseFoil];

figure('Color',[1 1 1]);
hold on
for i = 1:2
    switch i
        case {1}
            color = 'g';
        case {2}
            color = 'b';
    end 
shadedErrorBar([],fanceyFigMeans(i,:),fanceyFigStds(i,:),color,1) 
p(i) = shadedErrorBar([],fanceyFigMeans(i,:),fanceyFigStds(i,:),color,1);
end
title('Mean Pupil Size')
legend([p(1,1).mainLine p(1,2).mainLine],'Match','Foil');
xlabel('Time (ms)')
ylabel('Pupil Size (Normalized)')
vline(mean(probeStartms), 'r', 'Probe')
vline(mean(cueStartms), 'r', 'Cue')
vline(mean(delayStartms), 'r', 'Delay')
hold off

%% Plotting Pupil Size vs. viewing time during probe
% scatter plot and regression - viewing time of match (and foil and
% meither) vs. pupil size
r = corrcoef(meanPSDevProbe, propMatchTime');
figure;
subplot(2,2,1), scatter(meanPSDevProbe, propMatchTime'), lsline 
title(sprintf('Viewing time of match as a function of pupil size, r = %0.5g',r(2)))
ylabel('Match Viewing Time (s)')
xlabel('Changes in Pupil Size (A.U.)')
r = corrcoef(meanPSDevProbe, relativeMatchTime');
subplot(2,2,2), scatter(meanPSDevProbe, relativeMatchTime'), lsline 
title(sprintf('Relative Viewing time of match as a function of pupil size, r = %0.5g',r(2)))
ylabel('Relative viewing time: Match vs. Foil')
xlabel('Changes in Pupil Size (A.U.)')
r = corrcoef(meanPSAbsDevProbe, propMatchTime');
subplot(2,2,3), scatter(meanPSAbsDevProbe, propMatchTime'), lsline 
title(sprintf('Viewing time of match as a function of absolute change in pupil size, r = %0.5g',r(2)))
ylabel('Match Viewing Time (s)')
xlabel('Absolute Changes in Pupil Size (A.U.)')
r = corrcoef(meanPSAbsDevProbe, relativeMatchTime');
subplot(2,2,4), scatter(meanPSAbsDevProbe, relativeMatchTime'), lsline
title(sprintf('Relative viewing time of match as a function of absolute change in pupil size, r = %0.5g',r(2)))
ylabel('Relative viewing time: Match vs. Foil')
xlabel('Absolute Changes in Pupil Size (A.U.)')


% WITH Z-SCORES
r = corrcoef(meanPSZprobe, propMatchTime');
figure;
subplot(1,2,1), scatter(meanPSZprobe, propMatchTime'), lsline 
title(sprintf('Viewing time vs. pupil size during the probe, r = %0.5g',r(2)))
ylabel('Match Viewing Time (s)')
xlabel('Changes in Pupil Size (Z-Score)')
r = corrcoef(meanPSZprobe, relativeMatchTime');
subplot(1,2,2), scatter(meanPSZprobe, relativeMatchTime'), lsline 
title(sprintf('Relative Viewing time vs. pupil size during the probe, r = %0.5g',r(2)))
ylabel('Relative viewing time: Match vs. Foil')
xlabel('Changes in Pupil Size (Z-Score)')


%% logistic regression pupil size vs. choosing match vs. foil

figure;
[b,dev,stats] = glmfit(meanPSDevProbe,(relativeMatchTime > .5),'binomial','link','logit');
xx = linspace(min(meanPSDevProbe)+min(meanPSDevProbe)*.1,max(meanPSDevProbe)+max(meanPSDevProbe)*.1);
yfit = glmval(b,xx,'logit');
subplot(121), plot(meanPSDevProbe,(relativeMatchTime > .5),'o',xx,yfit,'-')
title('Probability of viewing Match or Foil more as a function of changes in pupil size')
ylabel('Probability that the subject will look at the Match more than the Foil')
xlabel('Changes in Pupil Size (A.U.)')

% absolute changes in pupil size 
[b,dev,stats] = glmfit(meanPSAbsDevProbe,(relativeMatchTime > .5),'binomial','link','logit');
xx = linspace(min(meanPSAbsDevProbe)+min(meanPSAbsDevProbe)*.1,max(meanPSAbsDevProbe)+max(meanPSAbsDevProbe)*.1);
yfit = glmval(b,xx,'logit');
subplot(122), plot(meanPSAbsDevProbe,(relativeMatchTime > .5),'o',xx,yfit,'-')
title('Probability of viewing Match or Foil more as a function of changes in pupil size')
ylabel('Probability that the subject will look at the Match more than the Foil')
xlabel('Absolute Changes in Pupil Size (A.U.)')

% With Z-score
figure;
[b,dev,stats] = glmfit(meanPSZprobe,(relativeMatchTime > .5),'binomial','link','logit');
xx = linspace(min(meanPSZprobe)+min(meanPSZprobe)*.1,max(meanPSZprobe)+max(meanPSZprobe)*.1);
yfit = glmval(b,xx,'logit');
plot(meanPSZprobe,(relativeMatchTime > .5),'o',xx,yfit,'-')
title('Probability of viewing Match or Foil more as a function of changes in pupil size during the Probe')
ylabel('Probability that the subject will look at the Match more than the Foil')
xlabel('Changes in Pupil Size (Z-Score)')


%% Boxplots of pupil size changes for match and foil during the probe

figure;
subplot(121), boxplot(meanPSDevProbe,relativeMatchTime > .5)
title(sprintf('Pupil Size - %s', dv.subj))
ylabel('Changes in pupil size (A.U.)')
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

subplot(122), boxplot(meanPSAbsDevProbe,relativeMatchTime > .5)
title(sprintf('Pupil Size - %s', dv.subj))
ylabel('Absolute Changes in pupil size (A.U.)')
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

% With Z-score
figure;
boxplot(meanPSZprobe,relativeMatchTime > .5)
title(sprintf('Pupil Size - %s', dv.subj))
ylabel('Changes in pupil size (Z-score)')
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);


%% logistic regression pupil size DURING SCENE CUE vs. choosing match vs. foil

figure;
[b,dev,stats] = glmfit(meanPSDevCue,(relativeMatchTime > .5),'binomial','link','logit');
xx = linspace(min(meanPSDevCue)+min(meanPSDevCue)*.1,max(meanPSDevCue)+max(meanPSDevCue)*.1);
yfit = glmval(b,xx,'logit');
subplot(121), plot(meanPSDevCue,(relativeMatchTime > .5),'o',xx,yfit,'-')
title('During Scene Cue')
ylabel('Probability that the subject will look at the Match more than the Foil')
xlabel('Changes in Pupil Size (A.U.)')

% absolute changes in pupil size 
[b,dev,stats] = glmfit(meanPSAbsDevCue,(relativeMatchTime > .5),'binomial','link','logit');
xx = linspace(min(meanPSAbsDevCue)+min(meanPSAbsDevCue)*.1,max(meanPSAbsDevCue)+max(meanPSAbsDevCue)*.1);
yfit = glmval(b,xx,'logit');
subplot(122), plot(meanPSAbsDevCue,(relativeMatchTime > .5),'o',xx,yfit,'-')
title('During Scene Cue')
ylabel('Probability that the subject will look at the Match more than the Foil')
xlabel('Absolute Changes in Pupil Size (A.U.)')

% Z-score
figure;
[b,dev,stats] = glmfit(meanPSZcue,(relativeMatchTime > .5),'binomial','link','logit');
xx = linspace(min(meanPSZcue)+min(meanPSZcue)*.1,max(meanPSZcue)+max(meanPSZcue)*.1);
yfit = glmval(b,xx,'logit');
plot(meanPSZcue,(relativeMatchTime > .5),'o',xx,yfit,'-')
title('During Scene Cue')
ylabel('Probability that the subject will look at the Match more than the Foil')
xlabel('Changes in Pupil Size (Z-Score)')


%% Boxplots of pupil size changes for match and foil during the SCENE CUE

figure;
subplot(121), boxplot(meanPSDevCue,relativeMatchTime > .5)
title(sprintf('Pupil Size during Scene Cue- %s', dv.subj))
ylabel('Changes in pupil size (A.U.)')
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

subplot(122), boxplot(meanPSAbsDevCue,relativeMatchTime > .5)
title(sprintf('Pupil Size during Scene Cue - %s', dv.subj))
ylabel('Absolute Changes in pupil size (A.U.)')
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);


% ANOVA
[P,ANOVATAB,STATS] = anova1(meanPSAbsDevCue, relativeMatchTime > .5);
title(sprintf('Reaction Time of First Saccade - Subject: %s, F = %0.5g', dv.subj,ANOVATAB{2,5}))
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
ylabel('Reaction Time (s)')

% Z-score
figure;
boxplot(meanPSZcue,relativeMatchTime > .5)
title(sprintf('Pupil Size during Scene Cue- %s', dv.subj))
ylabel('Changes in pupil size (Z-score)')
xtix = {'Match','Foil'};
xtixloc = [1 2];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);


%% Scatter plots - viewing vs pupil size during SCENE CUE

r = corrcoef(meanPSDevCue, propMatchTime');
figure;
subplot(2,2,1), scatter(meanPSDevCue, propMatchTime'), lsline 
title(sprintf('Viewing time of match as a function of pupil size during scene cue, r = %0.5g',r(2)))
ylabel('Match Viewing Time (s)')
xlabel('Changes in Pupil Size (A.U.)')
r = corrcoef(meanPSDevCue, relativeMatchTime');
subplot(2,2,2), scatter(meanPSDevCue, relativeMatchTime'), lsline 
title(sprintf('Relative Viewing time of match as a function of pupil size during scene cue, r = %0.5g',r(2)))
ylabel('Relative viewing time: Match vs. Foil')
xlabel('Changes in Pupil Size (A.U.)')
r = corrcoef(meanPSAbsDevCue, propMatchTime');
subplot(2,2,3), scatter(meanPSAbsDevCue, propMatchTime'), lsline 
title(sprintf('Viewing time of match as a function of absolute change in pupil size during scene cue, r = %0.5g',r(2)))
ylabel('Match Viewing Time (s)')
xlabel('Absolute Changes in Pupil Size (A.U.)')
r = corrcoef(meanPSAbsDevCue, relativeMatchTime');
subplot(2,2,4), scatter(meanPSAbsDevCue, relativeMatchTime'), lsline
title(sprintf('Relative viewing time of match as a function of absolute change in pupil size during scene cue, r = %0.5g',r(2)))
ylabel('Relative viewing time: Match vs. Foil')
xlabel('Absolute Changes in Pupil Size (A.U.)')

% Z-score
r = corrcoef(meanPSZcue, propMatchTime');
figure;
subplot(1,2,1), scatter(meanPSZcue, propMatchTime'), lsline 
title(sprintf('Viewing time of match as a function of pupil size during scene cue, r = %0.5g',r(2)))
ylabel('Match Viewing Time (s)')
xlabel('Changes in Pupil Size (Z-score)')
r = corrcoef(meanPSZcue, relativeMatchTime');
subplot(1,2,2), scatter(meanPSZcue, relativeMatchTime'), lsline 
title(sprintf('Relative Viewing time of match as a function of pupil size during scene cue, r = %0.5g',r(2)))
ylabel('Relative viewing time: Match vs. Foil')
xlabel('Changes in Pupil Size (Z-score)')


%% ploting pupil size during a particular trial
trialnum = 40;

%plot(pcSigChang{2})

figure; 
subplot(131),plot(pupilSizeDev{trialnum})
title('Deviations in Pupil Size (Mean = 0)')
xlabel('Time (ms)')
ylabel('A.U.')
subplot(132),plot(pupilSizeAbsDev{trialnum})
title('Absolute deviations in Pupil Size (Mean = 0)')
xlabel('Time (ms)')
ylabel('A.U.')
subplot(133),plot(pupilSizeZ{trialnum})
title('Deviations in Pupil Size')
xlabel('Time (ms)')
ylabel('Changes in Pupil Size (Z-score)')

%% Plotting during probe state during a particular trial
trialnum = 20;
figure; 
subplot(121),plot(psDevProbe{trialnum})
title('Deviations in Pupil Size (Mean = 0)')
xlabel('Time (ms)')
ylabel('A.U.')
subplot(122),plot(psAbsDevProbe{trialnum})
title('Absolute deviations in Pupil Size (Mean = 0)')
xlabel('Time (ms)')
ylabel('A.U.')


%% plot all trials for pupil size 
figure; 
n = 5; %length(pupilSizeDev);
for i = 1:n
subplot(1,n,i),plot(pupilSizeDev{i})
title('Deviations in Pupil Size (Mean = 0)')
xlabel('Time (ms)')
ylabel('A.U.')
end


%% Signal Detection Theory Analysis

%% Heatmaps overlayed on scenes


%% Save

%publish func or export_fig (needs ghostscript)
% add date
% export_fig dv.subj Summary -pdf -append
% 
% publish(dv.subj, 'pdf')






