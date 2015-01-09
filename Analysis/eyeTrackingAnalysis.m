% Analyze eye tracking data
% PDS.eyepos = [x,y,trialTime, trialState], each cell from a trial
% PDS.data.eyeposloop is the same i think
% dv.states.SHOWPROBE = 11;
% objLocs = % [correctObj, foilObj, loc1, loc2, correctLoc]
% PDS.data.correctObject(iTrial)
   
    format longg
    
    eyePos = PDS.eyepos(strcmp(PDS.trialType,'test'));
    objLocs = PDS.data.objectLocs(strcmp(PDS.trialType,'test'));
    
    state = dv.states.SHOWPROBE; %only the eyePos during memory probe (just in case we want to switch what state we are analyzing)
    
    % NOTE: Double check that the foil obj is randomly paired with each stimulus each trial
    
    % NOTE: zeros in cells in eyePosProp are from break fixation, there are no eye positions during the show probe state then
    
    n = length(eyePos);
    eyePosState = cell(1,n);
    eyePosProp = cell(1,n);
    emptyCells = zeros(n,1);
    firstSaccadeSample = zeros(n,1);
    firstSaccadeTimeInState = zeros(n,1);
    firstSaccadeStatus = zeros(n,1);
    relativeMatchTime = zeros(n,1);
    propMatchTime = zeros(n,1);
    propFoilTime = zeros(n,1);
    propNeitherTime = zeros(n,1);
    
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
            eyePosStateTime = eyePosState{iTrial}(:,3);
            firstSaccadeTimeInState(iTrial) = eyePosStateTime(firstSaccadeSample(iTrial)) - eyePosStateTime(1);
            % status of first saccade (location - target or foil)
            firstSaccadeStatus(iTrial) = eyePosProp{iTrial}(firstSaccadeSample(iTrial));
            
            %%%%%%%%%%%%%%%%%%%%%%
            
            % relative viewing
            sampleTime = (eyePosStateTime(end) - eyePosStateTime(1)) / length(eyePosStateTime);
            propMatchTime(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 1)); % why are these all the same??
            propFoilTime(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 2));
            propNeitherTime(iTrial) = sampleTime * length(find(eyePosProp{iTrial} == 0));
            relativeMatchTime(iTrial) = propMatchTime(iTrial) / (propMatchTime(iTrial) + propFoilTime(iTrial));
        else
            firstSaccadeSample(iTrial) = NaN;
            firstSaccadeTimeInState(iTrial) = NaN;
            firstSaccadeStatus(iTrial) = NaN;
            propMatchTime(iTrial) = NaN;
            propFoilTime(iTrial) = NaN;
            propNeitherTime(iTrial) = NaN;
            relativeMatchTime(iTrial) = NaN;
        end
        
        
    end
    
    % eliminate NaNs
    firstSaccadeSample = firstSaccadeSample(~isnan(firstSaccadeSample));
    firstSaccadeTimeInState = firstSaccadeTimeInState(~isnan(firstSaccadeTimeInState));
    firstSaccadeStatus = firstSaccadeStatus(~isnan(firstSaccadeStatus));
    propMatchTime = propMatchTime(~isnan(propMatchTime));
    propFoilTime = propFoilTime(~isnan(propFoilTime));
    propNeitherTime = propNeitherTime(~isnan(propNeitherTime));
    relativeMatchTime = relativeMatchTime(~isnan(relativeMatchTime));
    

%% Plotting
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

%%
% RTs
saccadeRT = [firstSaccadeTimeInState firstSaccadeStatus];
anova1(saccadeRT(:,1),saccadeRT(:,2))


firstSaccadeTimeInState = sort(firstSaccadeTimeInState);  % sorting we lose what the 'status' is
firstSaccadeTimeInState1 = firstSaccadeTimeInState(1:end/2);
firstSaccadeTimeInState2 = firstSaccadeTimeInState(end/2+1:end);


%% logistic regression on RTs and first saccades
saccadeRT = [firstSaccadeTimeInState firstSaccadeStatus];
x = saccadeRT(:,1);
firstSaccadeStatus(firstSaccadeStatus == 2) = 0; % need a binary variable
y = firstSaccadeStatus;

figure;
[b,dev,stats] = glmfit(x,y,'binomial','link','logit');
yfit = glmval(b,x,'logit');

xx = linspace(0,1.8);
yfit = glmval(b,xx,'logit');
plot(x,y,'o',xx,yfit,'-')









