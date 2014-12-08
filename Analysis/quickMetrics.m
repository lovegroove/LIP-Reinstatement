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

%%
% oddsCorrectSession(125:164)
% mean(oddsCorrectSession(125:164))
propNonmatchTime = dv.pa.probeTime - (propMatchTime(125:164) + prop0Session(125:164));
y = [propMatchTime(125:164) prop0Session(125:164) propNonmatchTime]; %specific indices here NOTE
hold on
bar(mean(y))
errorbar(mean(y),std(y), 'r', 'Marker', 'none', 'LineStyle', 'none' );
hold off

hold on
subplot(121), bar(mean(y)), errorbar(mean(y),std(y),'.')
subplot(122), boxplot(y)
hold off

subplot(1,2,1), boxplot(propMatchSession(125:164)) % maybe call it prop matching instead
subplot(1,2,2), hist(propMatchSession(125:164))
mean(propMatchSession(125:164))
std(propMatchSession(125:164))
median(propMatchSession(125:164))

propTrialsCorrect = length(find(propMatchSession(125:164) > .5)) / length(propMatchSession(125:164)); 
length(find(propMatchSession(125:164) < .5)) 

firstSaccadeCorrect(125:164)

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