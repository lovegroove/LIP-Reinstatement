% Analyze eye tracking data
% PDS.eyepos = [x,y,trialTime, trialState], each cell from a trial
% PDS.data.eyeposloop is the same i think
% dv.states.SHOWPROBE = 11;

format longg

eyePos = PDS.eyepos(strcmp(PDS.trialType,'test'));
objLocs = PDS.data.objectLocs(strcmp(PDS.trialType,'test'));

state = dv.states.SHOWPROBE; %only the eyePos during memory probe

% NOTE: Double check that the objlocs always have the correct target first in the array!!!!!!!!!!!!!!!!!!!! 

n = length(eyePos);
eyePosState = cell(1,n);
eyePosProp = cell(1,n);
emptyCells = zeros(n,1);
firstSaccadeSample = zeros(n,1);
firstSaccadeTimeInState = zeros(n,1);
firstSaccadeStatus = zeros(n,1);
relativeMatch = zeros(n,1);
relativeMatchTime = zeros(n,1);

% save these as awe go???
prop1 = zeros(n,1);
prop2 = zeros(n,1);
prop0 = zeros(n,1);
prop1Time = zeros(n,1);
prop2Time = zeros(n,1);
prop0Time = zeros(n,1);

for iTrial = 1:n
    
    eyePosState{iTrial} = eyePos{iTrial}(eyePos{iTrial}(:,4)==state,:); % only eye position during a certain state
    sampleLen = length(eyePosState{iTrial});
    
    for jSample = 1:sampleLen
        % eyes within correct object boundaries
        if eyePosState{iTrial}(jSample,1) >= objLocs{iTrial}{2,1}(1,1) && eyePosState{iTrial}(jSample,1) <= objLocs{iTrial}{2,1}(1,3) && eyePosState{iTrial}(jSample,2) >= objLocs{iTrial}{2,1}(1,2) && eyePosState{iTrial}(jSample,2) <= objLocs{iTrial}{2,1}(1,4)
            
        eyePosProp{iTrial}(jSample,1) = 1;    
            
        % eyes within foil object boundaries
        elseif eyePosState{iTrial}(jSample,1) >= objLocs{iTrial}{2,2}(1,1) && eyePosState{iTrial}(jSample,1) <= objLocs{iTrial}{2,2}(1,3) && eyePosState{iTrial}(jSample,2) >= objLocs{iTrial}{2,2}(1,2) && eyePosState{iTrial}(jSample,2) <= objLocs{iTrial}{2,2}(1,4)
            
        eyePosProp{iTrial}(jSample,1) = 2;  
        
        % eyes in neither
        else
        eyePosProp{iTrial}(jSample,1) = 0;  
        
        end
    end
    
    % flag empty cells due to break fixation
    emptyCells(iTrial) = isempty(eyePosProp{iTrial});
    
    
    if ~isempty(eyePosProp{iTrial})
        % first saccade
        firstSaccadeSample(iTrial) = find(eyePosProp{iTrial} ~= 0,1);
        % time of first saccade
        eyePosStateTime = eyePosState{iTrial}(:,3);
        firstSaccadeTimeInState(iTrial) = eyePosStateTime(firstSaccadeSample(iTrial)) - eyePosStateTime(1);
        % status of first saccade (location - target or foil)
        firstSaccadeStatus(iTrial) = eyePosProp{iTrial}(firstSaccadeSample(iTrial)); 
        
        % relative viewing
        prop1(iTrial) = length(find(eyePosProp{iTrial} == 1)) / length(eyePosProp{iTrial});
        prop2(iTrial) = length(find(eyePosProp{iTrial} == 2)) / length(eyePosProp{iTrial});
        prop0(iTrial) = length(find(eyePosProp{iTrial} == 0)) / length(eyePosProp{iTrial});
        
        sampleTime = (eyePosStateTime(end) - eyePosStateTime(1)) / length(eyePosStateTime);
        prop1Time(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 1));
        prop2Time(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 2));
        prop0Time(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 0));
        
        relativeMatch(iTrial) = prop1(iTrial) / (prop1(iTrial) + prop2(iTrial));
        relativeMatchTime(iTrial) = prop1Time(iTrial) / (prop1Time(iTrial) + prop2Time(iTrial));
    else
        firstSaccadeSample(iTrial) = NaN;
        firstSaccadeTimeInState(iTrial) = NaN;
        firstSaccadeStatus(iTrial) = NaN;
    end
    
    
end

%eyePosProp = eyePosProp(~cellfun('isempty',eyePosProp));  

% eliminate NaNs
firstSaccadeSample = firstSaccadeSample(~isnan(firstSaccadeSample));
firstSaccadeTimeInState = firstSaccadeTimeInState(~isnan(firstSaccadeTimeInState));
firstSaccadeStatus = firstSaccadeStatus(~isnan(firstSaccadeStatus));




%% Analysis


%% Plotting
% percent correct
pcTrials = length(find(relativeMatch > .5)) / length(relativeMatch);

figure;
% Total viewing time in a session
relativeNonmatchTime = dv.pa.probeTime - (relativeMatchTime + prop0Time); %?????? can't have negatives
y = [relativeMatchTime relativeNonmatchTime prop0Time];
subplot(221), h = barwitherr(std(y),mean(y));
title(sprintf('Total Viewing time - %s', dv.subj)) 
set(gca,'XTickLabel',{'Matching','Non-matching','Neither'})
set(h,'FaceColor','b');
ylabel('Viewing Time (s)')
subplot(222), boxplot(y)
title(sprintf('Total Viewing time - %s', dv.subj))
ylabel('Viewing Time (s)')
xtix = {'Matching','Non-matching','Neither'};
xtixloc = [1 2 3];
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

%  Relative Viewing matched vs. non-matched
subplot(2,2,3), boxplot(relativeMatch)
title(sprintf('Proportion of Viewing Time: matching vs. the non-matching, Percent Correct: %0.5g',pcTrials * 100))
ylabel('Relative Viewing Time: matching vs. the non-matching')
subplot(2,2,4), hist(relativeMatch)
g = findobj(gca,'Type','patch');
set(g,'FaceColor','b');
title(sprintf('Relative Viewing Time: matching vs. the non-matching, Percent Correct: %0.5g',pcTrials * 100))
xlabel('Relative Viewing Time: matching vs. the non-matching')
ylabel('Number of Trials')




%%
% RTs
saccadeRT = [firstSaccadeTimeInState firstSaccadeStatus];
anova1(saccadeRT(:,1),saccadeRT(:,2))


firstSaccadeTimeInState = sort(firstSaccadeTimeInState);  % sorting we lose what the 'status' is
firstSaccadeTimeInState1 = firstSaccadeTimeInState(1:end/2);
firstSaccadeTimeInState2 = firstSaccadeTimeInState(end/2+1:end);


%% logistic regression on RTs and first saccades
x = saccadeRT(:,1);
firstSaccadeStatus(firstSaccadeStatus == 2) = 0; % need a binary variable
y = firstSaccadeStatus;

figure;
[b,dev,stats] = glmfit(x,y,'binomial','link','logit');
yfit = glmval(b,x,'logit');

xx = linspace(0,1.8);
yfit = glmval(b,xx,'logit');
plot(x,y,'o',xx,yfit,'-')









