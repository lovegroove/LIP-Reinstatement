function dv = liprein(dv)
% Condition file for LIP Reinstatement experiemnet, linked to
% runLIPReinstatement, and called by runPLDAPScv (cv = color version, which
% is needed)

% Trial Function & Trial Type
dv.trialFunction = 'runLIPreinTrial';
dv.filePaths = 'shopRig'; % set here rather than in loadImages func
dv.trialType = 'study'; % choose trial type, 'study' or 'test'
dv.finish = 1e3; % # of trials
dv.singleSession = 1;
dv.pa.singleSessionStudy = 40;
dv.pa.singleSessionTest = 80;
dv.pa.strictDelay = 1; %boolean

dv = pdsDefaultTrialStructure(dv); % calls dv = defaultColors(dv) IMPORTANT for assigning CLUT values; current overwriting some defaults too

Datapixx('Open') % HACK: Must reopen datapixx beacause I'm not using overlay so I can get colors

dv.disp.preflipbuffer = 10e-3; % 10 ms preflip (lots of textures to draw), need?

%-------------------------------------------------------------------------%
% Options - both of these must be flagged if you want to skip eyetracker
% completely and have it run right
dv.pass = 0; % ignore eyetracker by flagging fixationHeld as always true

dv.mouse = 0; % just toggle this if you want to use the mouse
if dv.mouse % this is dumb i know, just to quickly simplify
    dv.useMouse = 1; % sets cursor to eyes via pdsGetEyePosition func
    dv.showMouse = 1;
else
    dv.useMouse = 0; % sets cursor to eyes via pdsGetEyePosition func
    dv.showMouse = 0;
end
    
%-------------------------------------------------------------------------%
% Don't forget to CHANGE FILE PATHS on other rigs
% Load images function - currently calls new pairs on every "new" session
% So remember that!
%****************** 
if dv.newsession 
dv = loadImages(dv);
end
%*******************
%-------------------------------------------------------------------------%
%% Parameters - dv.pa

% Stimulus Size and Location (dv.disp - display parameters)
dv.pa.objectSize = 4; % in degrees

% step size around the center, angle at which stimulus is presented
dv.pa.stimThetas = 0:45:360;

%%%%% Test Trial Stimulus Location and Geometry Parameters
dv.pa.alpha = 8; % degrees of visual angle
[dv.pa.Dx, dv.pa.Dy] = calcVisAngDS(dv.pa.alpha, dv.disp.viewdist, dv.disp.widthcm, dv.disp.heightcm, dv.disp.winRect(3), dv.disp.winRect(4)); % Func to calculate visual angle

% Get the center coordinate of the window
[dv.pa.xCenter, dv.pa.yCenter] = RectCenter(dv.disp.winRect);

% Make Delay box (no longer using delay box)
delayRect = [0 0 50 50];
dv.pa.centeredDelayRect = CenterRectOnPointd(delayRect, dv.pa.xCenter, dv.pa.yCenter);
dv.pa.delayBoxColor = [0, 0, 0]; % delay box color

%% Time & Space

% Task Events - (was in ms, keep it that way?)
dv.pa.sceneTime = 1;
dv.pa.showPairTime = 2;
dv.pa.probeCueTime = .75;
dv.pa.delayTime = 4.25; % delay and probe cue time should add up to your desired total delay
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
xCoords = [-dv.pa.fixCrossDimPix dv.pa.fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero, center set automatically)
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
dv.states.TRIALCOMPLETE = 6;
dv.states.BREAKFIX  = 7;

dv.states.FPONDELAY = 12;
dv.states.FPHOLDDELAY = 13;

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


