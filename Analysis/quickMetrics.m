% Quick simple analysis of eye tracking data to see if we are observing a
% behavioral effect
%

n = length(PDS.data.eyeLocProbe);
oddsMatchSession = zeros(n,1);
propMatchSession = zeros(n,1);
prop1Session = zeros(n,1);
prop2Session = zeros(n,1);
prop0Session = zeros(n,1);
firstSaccade = zeros(n,1);
firstSaccadeCorrect = zeros(n,1);
propMatchTime = zeros(n,1);


for i = 1:n 
    loopTime = dv.pa.probeTime / length(PDS.data.eyeLocProbe{i});
    
    prop1 = length(find(PDS.data.eyeLocProbe{i} == 1)) / length(PDS.data.eyeLocProbe{i});
    prop2 = length(find(PDS.data.eyeLocProbe{i} == 2)) / length(PDS.data.eyeLocProbe{i});
    prop0 = length(find(PDS.data.eyeLocProbe{i} == 0)) / length(PDS.data.eyeLocProbe{i});
    
    prop1Time = loopTime * length(find(PDS.data.eyeLocProbe{i} == 1));
    prop2Time = loopTime * length(find(PDS.data.eyeLocProbe{i} == 2));
    prop0Time = loopTime * length(find(PDS.data.eyeLocProbe{i} == 0));
    
    if PDS.data.correctObject(i) == 1
        %oddsMatch = prop1 / (1-prop1); % not including looking at neither
        propMatch = prop1 / (prop1 + prop2); % shouldn't this be the other way
        propMatchTime(i) = prop1Time;
    elseif PDS.data.correctObject(i) == 2
        %oddsMatch = prop2 / (1-prop2); 
        propMatch = prop2 / (prop1 + prop2);
        propMatchTime(i) = prop2Time;
        
    else
        propMatch = NaN;
        propMatchTime(i) = NaN;
    end
    
    %oddsMatchSession(i) = oddsMatch;
    propMatchSession(i) = propMatch;
    prop1Session(i) = prop1;
    prop2Session(i) = prop2;
    prop0Session(i) = prop0;
    
    % first saccades
    if find(PDS.data.eyeLocProbe{i} == 1,1) < find(PDS.data.eyeLocProbe{i} == 2,1)
        firstSaccade(i) = 1;
        if PDS.data.correctObject(i) == 1 && firstSaccade(i) == 1
            firstSaccadeCorrect(i) = 1;
        else
            firstSaccadeCorrect(i) = 0;
        end
    elseif find(PDS.data.eyeLocProbe{i} == 2,1) < find(PDS.data.eyeLocProbe{i} == 1,1)
        firstSaccade(i) = 2;
        if PDS.data.correctObject(i) == 2 && firstSaccade(i) == 2
            firstSaccadeCorrect(i) = 1;
        else
            firstSaccadeCorrect(i) = 0;
        end
    end

end

%%  Plotting

trialRange = [123 162];

% Total viewing time in a session
propNonmatchTime = dv.pa.probeTime - (propMatchTime(trialRange(1):trialRange(2)) + prop0Session(trialRange(1):trialRange(2))); %125:164
y = [propMatchTime(trialRange(1):trialRange(2)) propNonmatchTime prop0Session(trialRange(1):trialRange(2))]; %specific indices here NOTE , propNon is alreayd specific indicies
subplot(221), h = barwitherr(std(y),mean(y));
title('Total Viewing time - 1 Session')
set(gca,'XTickLabel',{'Matching','Non-matching','Neither'})
set(h,'FaceColor','b');
ylabel('Viewing Time (s)')
subplot(222), boxplot(y)
title('Total Viewing time - 1 Session')
ylabel('Viewing Time (s)')
xtix = {'Matching','Non-matching','Neither'};   
xtixloc = [1 2 3];      
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

%  Proportion viewing matched vs. non-matched
subplot(2,2,3), boxplot(propMatchSession(trialRange(1):trialRange(2))) 
title('Proportion of time viewing the matched vs. the non-matched stimulus')
ylabel('Proportion of time viewing the matching stimulus vs. non-matching')
subplot(2,2,4), hist(propMatchSession(trialRange(1):trialRange(2)))
g = findobj(gca,'Type','patch');
set(g,'FaceColor','b');
title('Proportion of time viewing the matched vs. the non-matched stimulus')
xlabel('Proportion of time viewing the matching stimulus vs. non-matching')
ylabel('Number of Trials')
% mean(propMatchSession(125:164))
% std(propMatchSession(125:164))
% median(propMatchSession(125:164))

propTrialsCorrect = length(find(propMatchSession(trialRange(1):trialRange(2)) > .5)) / length(propMatchSession(trialRange(1):trialRange(2)));
propTrialsCorrect

%% Across Sessions

for k = 1:length(propMatchTime)
    if isnan(propMatchTime(k))
        propMatchTime(k) = [];
    end
    if isnan(prop0Session(k))
        prop0Session(k) = [];
    end
    
end

% Total viewing time 
propNonmatchTime = dv.pa.probeTime - (propMatchTime + prop0Session); 
y = [propMatchTime propNonmatchTime prop0Session]; 
subplot(221), h = barwitherr(std(y),mean(y));
title('Total Viewing time - Across Sessions')
set(gca,'XTickLabel',{'Matching','Non-matching','Neither'})
set(h,'FaceColor','b');
ylabel('Viewing Time (s)')
subplot(222), boxplot(y)
title('Total Viewing time - Across Sessions')
ylabel('Viewing Time (s)')
xtix = {'Matching','Non-matching','Neither'};   
xtixloc = [1 2 3];      
set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

subplot(2,2,3), boxplot(propMatchSession) 
title('Proportion of time viewing the matched vs. the non-matched stimulus')
ylabel('Proportion of time viewing the matching stimulus vs. non-matching')
subplot(2,2,4), hist(propMatchSession);
g = findobj(gca,'Type','patch');
set(g,'FaceColor','b');
title('Proportion of time viewing the matched vs. the non-matched stimulus')
xlabel('Proportion of time viewing the matching stimulus vs. non-matching')
ylabel('Number of Trials')

%%

% plot(propMatchSession(125:164))
% 
prop1Session(125:164)
prop2Session(125:164)
prop0Session(125:164)

x = [ones(length(propMatchSession(125:164)),1) firstSaccadeCorrect(125:164)];
[B,BINT,R,RINT,STATS] = regress(propMatchSession(125:164),x);
regstats(propMatchSession(125:164),firstSaccadeCorrect(125:164))
scatter(yhat,r)

%---------------------
% single trial
% proportion of viewing time over trial
% heat map superimposed on images