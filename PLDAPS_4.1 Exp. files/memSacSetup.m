function p=memSacSetup(p)
% Memory guided saccades experiment setup
% example:
% load settingsStruct
% p=pldaps(@memSacSetup,'subj', settingsStruct);
% p.run

% Trial function
p.defaultParameters.pldaps.trialFunction='memSac';

% Default trial parameters (stimulus stuff etc.)
p = pdsDefaultTrialStructure(p); 
    
%         dv.defaultParameters.pldaps.trialMasterFunction='runTrial'; % automatically specified, right?

% trial length (s)
p.trial.pldaps.maxTrialLength = 5;
p.trial.pldaps.maxFrames = p.trial.pldaps.maxTrialLength*p.trial.display.frate;

c.Nr=1; % number of conditions
p.conditions=repmat({c},1,200);

p.defaultParameters.pldaps.finish = length(p.conditions);

% Initializes all interated trial states
defaultTrialVariables(p);

% Reward
p.trial.stimulus.rewardTime = .1;

% Fixation
p.trial.stimulus.preTrial     = .5;       
p.trial.stimulus.fixWait      = 3;        
p.trial.stimulus.fixHold      = 1;
p.trial.stimulus.fixationXY   = [0 0]; % degrees

p.trial.stimulus.fpOffset = .5; % adding this in... unsure of value
%dv.trial.stimulus.fpOffset ??? needed?

% Targets
p.trial.stimulus.targWait   	= 1.5;
p.trial.stimulus.targHold   	= 0.5;
p.trial.stimulus.targOnset  	= [0.1 0.1];
p.trial.stimulus.targDuration 	= [2 .2];
p.trial.stimulus.targ1Loc        = [10 0]; % degrees (elh)
p.trial.stimulus.targ2Loc        = [-10 0]; % degrees (elh)

% Task
p.trial.stimulus.targUser  = 0;
p.trial.stimulus.proportionMemory = .5; % generates rate at which memory saccades are shown

p.trial.stimulus.targRange = [-20 -10 20 10]; 
p.trial.stimulus.jitterspace = 2;

% Windows
p.trial.stimulus.fpWin 			= [.8 .8]; % x,y radius
p.trial.stimulus.winScaleVisual = .2; % scale targwin with eccentricity
p.trial.stimulus.winScaleMemory = .4;

%-------------------------------------------------------------------------%
% Colors
p.defaultParameters.display.humanCLUT(12,:)=[1 0 0];
p.defaultParameters.display.monkeyCLUT(12,:)=[0 1 0];

p.defaultParameters.display.humanCLUT(13,:)=[0 1 0];
p.defaultParameters.display.monkeyCLUT(13,:)=p.defaultParameters.display.bgColor;
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.hGreen=12*[1 1 1]';
else
    p.defaultParameters.display.clut.hGreen=p.defaultParameters.display.humanCLUT(12+1,:)';
end

p.defaultParameters.display.humanCLUT(14,:)=[1 0 0];
p.defaultParameters.display.monkeyCLUT(14,:)=p.defaultParameters.display.bgColor;
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.hRed=13*[1 1 1]';
else
    p.defaultParameters.display.clut.hRed=p.defaultParameters.display.humanCLUT(13+1,:)';
end

p.defaultParameters.display.humanCLUT(15,:)=[0 0 0];
p.defaultParameters.display.monkeyCLUT(15,:)=p.defaultParameters.display.bgColor;
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.hBlack=14*[1 1 1]';
else
    p.defaultParameters.display.clut.hBlack=p.defaultParameters.display.humanCLUT(14+1,:)';
end

p.defaultParameters.display.humanCLUT(16,:)=[1 0 0];
p.defaultParameters.display.monkeyCLUT(16,:)=[1 0 0];
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.bRed=15*[1 1 1]';
else
    p.defaultParameters.display.clut.bRed=p.defaultParameters.display.humanCLUT(15+1,:)';
end

p.defaultParameters.display.humanCLUT(17,:)=[0 1 0];
p.defaultParameters.display.monkeyCLUT(17,:)=[0 1 0];
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.bGreen=16*[1 1 1]';
else
    p.defaultParameters.display.clut.bGreen=p.defaultParameters.display.humanCLUT(16+1,:)';
end

p.defaultParameters.display.humanCLUT(18,:)=[1 1 1];
p.defaultParameters.display.monkeyCLUT(18,:)=[1 1 1];
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.bWhite=17*[1 1 1]';
else
    p.defaultParameters.display.clut.bWhite=p.defaultParameters.display.humanCLUT(17+1,:)';
end

p.defaultParameters.display.humanCLUT(19,:)=[0 0 0];
p.defaultParameters.display.monkeyCLUT(19,:)=[0 0 0];
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.bBlack=18*[1 1 1]';
else
    p.defaultParameters.display.clut.bBlack=p.defaultParameters.display.humanCLUT(18+1,:)';
end

p.defaultParameters.display.humanCLUT(20,:)=[0 0 1];
p.defaultParameters.display.monkeyCLUT(20,:)=[0 0 1];
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.bBlue=19*[1 1 1]';
else
    p.defaultParameters.display.clut.bBlue=p.defaultParameters.display.humanCLUT(19+1,:)';
end

p.defaultParameters.display.humanCLUT(21,:)=[0 0 1];
p.defaultParameters.display.monkeyCLUT(21,:)=p.defaultParameters.display.bgColor;
if p.defaultParameters.datapixx.use && p.defaultParameters.display.useOverlay
    p.defaultParameters.display.clut.hBlue=20*[1 1 1]';
else
    p.defaultParameters.display.clut.hBlue=p.defaultParameters.display.humanCLUT(20+1,:)';
end


end