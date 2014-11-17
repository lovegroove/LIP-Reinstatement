function dv = liprein(dv)
% Condition file for LIP Reinstatement experiemnet, linked to
% runLIPReinstatement, and called by runPLDAPS

% Trial Function & Trial Type
dv.trialFunction = 'runLIPreinTrial';
dv.filePaths = 'shopRig'; % set here rather than in loadImages func now
dv.trialType = 'test'; % choose trial type, 'study' or 'test'
dv.finish = 1e3; % # of trials

dv = pdsDefaultTrialStructure(dv); % calls dv = defaultColors(dv) IMPORTANT for assigning CLUT values; also i am currently overwriting some dv.pa and dv.states below which could be cleaned up and customized

Datapixx('Open') % HACK: Must reopen datapixx beacause I'm not using overlay so I can get colors

dv.disp.preflipbuffer = 10e-3; % 10 ms preflip (lots of textures to draw), need????

% Options
dv.pass = 1; % ignore eyetracker by flagging fixationHeld as always true
dv.useMouse = 1; % sets cursor to eyes via pdsGetEyePosition func


% don't forget to CHANGE FILE PATHS
% Load images function - *******************right now, this is setup to call new pairings every time you run the experiment, so you have to be careful comparing study and test, because you might overwrite it with a new order
% need a better way to do this and to save it the the PDS struct in run
% file per trial...
dv = loadImages(dv);
% Proliferate pairings for many trials - 40 is the current size of a
% stimulus set: 8 shapes, each paired with 5 scenes
for i = 1:dv.finish / 40 
    dv.pairOrder= [dv.pairOrder; datasample(dv.pairOrder,40,1,'replace',false)];  %should we overwrite dv.pairOrder or no?
end

%% Parameters - dv.pa

% Stimulus Size and Location (dv.disp - display parameters)
dv.pa.objectSize = .3; % in percent of baseRect, better way to do this? 

%%%%% Test Trial Stimulus Location and Geometry Parameters
dv.pa.alpha = 15; % degrees of visual angle
[dv.pa.Dx, dv.pa.Dy] = calcVisAngDS(dv.pa.alpha, dv.disp.viewdist, dv.disp.widthcm, dv.disp.heightcm, dv.disp.winRect(3), dv.disp.winRect(4)); % Func to calculate visual angle

% Get the center coordinate of the window
[dv.pa.xCenter, dv.pa.yCenter] = RectCenter(dv.disp.winRect);

% Make Delay box
delayRect = [0 0 50 50];
dv.pa.centeredDelayRect = CenterRectOnPointd(delayRect, dv.pa.xCenter, dv.pa.yCenter);
dv.pa.delayBoxColor = [0, 0, 0]; % delay box color

%% Time & Space

% Task Events - (was in ms, keep it that way?)
dv.pa.sceneTime = 1;
dv.pa.showPairTime = 2;
dv.pa.probeCueTime = .5;
dv.pa.delayTime = 4.5; % delay and probe cue time should add up to your desired total delay
dv.pa.probeTime = 2;
dv.pa.graceTime = .5; % long enough?

% Reward time: is time that solenoid is opened for. set to 100 miliseconds
dv.pa.freeOn = 1; % if free trial, give reward independent of choice, right now we have no choice
dv.pa.rewardTime = .1;
dv.pa.rewardWait = 0; 
dv.pa.breakFixPenalty = 2;
dv.pa.jitterSize = .5;

% fixation
dv.pa.preTrial = .5;  
dv.pa.fixPtOffset = 1;
dv.pa.fixWait = 4;        
dv.pa.fixHold = 3; % 1
dv.pa.fixationXY = [0 0]; % degrees

%fixation cross
% Setup the text type for the window - surprised I can call screen here
% already
Screen('TextFont', dv.disp.ptr, 'Ariel');
Screen('TextSize', dv.disp.ptr, 36);

% Set the size of the arms of our fixation cross
dv.pa.fixCrossDimPix = 40;
xCoords = [-dv.pa.fixCrossDimPix dv.pa.fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero we will let the drawing routine center the cross in the center of our monitor for us)
yCoords = [0 0 -dv.pa.fixCrossDimPix dv.pa.fixCrossDimPix];
dv.pa.allCoords = [xCoords; yCoords];
dv.pa.lineWidthPix = 4; % line width for our fixation cross
dv.pa.fixCrossColor = [0, 0, 0]; % Fixation cross color


% window sizes 
dv.pa.fpWin = [1.5 1.5]; % x,y radius
dv.pa.winScale = .2; % scale targwin with eccentricity

% dot sizes for drawing
dv.pa.eyeW      = 8;    % eye indicator width in pixels
dv.pa.fixdotW   = 8;    % width of the fixation dot
dv.pa.targdotW  = 8;    % width of the target dot
dv.pa.cursorW   = 8;   % cursor width in pixels



%% States - dv.states (arbitrary values)

dv.states.START     = 1;
dv.states.FPON      = 2;
dv.states.FPHOLD    = 3;

%dv.states.VIRGINTRIAL = 4;

dv.states.TRIALCOMPLETE = 6;
dv.states.BREAKFIX  = 7;

dv.states.SHOWCUE = 8;
dv.states.SHOWPAIR = 9;
dv.states.DELAY = 10;
dv.states.SHOWPROBE = 11;


    %% PsychPort Audio
    dv = pdsAudioSetup(dv);
    
    %% Eyelink toolbox
    % Eyelink sampling rate, max trial length(s)
    dv = pdsEyelinkSetup(dv);
    

end


