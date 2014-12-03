%% Crude pass at plotting eyetracking data for LIP Reinstatement Experiment

% just playing

% Plot eye postion for each trial
figure;
hold on;
for k = 1:length(PDS.eyepos)-6 % currently just doing first 12
    
    xlim([0 dv.disp.winRect(3)])
    ylim([0 dv.disp.winRect(4)])
    for i = 1:length(PDS.eyepos{k})
        hold on;
        subplot(4,(length(PDS.eyepos)-6)/4,k), plot(PDS.eyepos{k}(i,1),PDS.eyepos{k}(i,2),'g--o','LineWidth',2)
    end
    
end
hold off;





%%
% One Trial - eye data
figure;
xlim([0 dv.disp.winRect(3)]) 
ylim([0 dv.disp.winRect(4)]) 
hold on
k = 2; % trial number
for i = 1:length(PDS.eyepos{k})
    plot(PDS.eyepos{k}(i,1),PDS.eyepos{k}(i,2),'g--o','LineWidth',2)
end
hold off

