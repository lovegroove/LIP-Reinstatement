function dv = liprein(dv)
% Condition file for LIP Reinstatement experiemnet, linked to
% runLIPReinstatement, and called by runPLDAPS

dv =[]; % normally in the runPLDAPS file but for testing purposes

% Trial Function & Trial Type
dv.trialFunction = 'runLIPreinTrial';
dv.trialType = 'study'; % choose trial type, 'study' or 'test'
dv.finish = 1e3; % # of trials
%dv = pdsDefaultTrialStructure(dv); % need/want this??????
dv.disp.preflipbuffer = 10e-3; % 10 ms preflip (lots of textures to draw), need????

% Options
dv.pass = 0; % ignore eyetracker
dv.useMouse = 0; % same as above? 


% don't forget to CHANGE FILE PATHS
% Load images function - *******************right now, this is setup to call new pairings every time you run the experiment, so you have to be careful comparing study and test, because you might overwrite it with a new order*********************************
% need a better way to do this and to save it the the PDS struct in run
% file per trial...
dv = loadImages(dv);
% Proliferate pairings for many trials - 40 is the current size of a
% stimulus set: 8 shapes, each paired with 5 scenes
for i = 1:dv.finish / 40 
    dv.pairOrder= [dv.pairOrder; datasample(dv.pairOrder,40,1,'replace',false)];  %should we overwrite dv.pairOrder or no?
end

%% Parameters - dv.pa

% Stimulus Size and Location (ds - display parameters)
dv.pa.objectSize = .25; % in percent of baseRect, better way to do this? 

%%%%% Test Trial Stimulus Location and Geometry Parameters
dv.pa.alpha = 6; % degrees of visual angle
[dv.pa.Dx, dv.pa.Dy] = calcVisAngDS(dv.pa.alpha, ds.viewdist, ds.widthcm, ds.heightcm, ds.windowRect(3), ds.windowRect(4)); % Func to calculate visual angle, relies on ds

% Get the center coordinate of the window
[dv.pa.xCenter, dv.pa.yCenter] = RectCenter(ds.windowRect);

% Make Delay box
delayRect = [0 0 50 50];
dv.pa.centeredDelayRect = CenterRectOnPointd(delayRect, dv.pa.xCenter, dv.pa.yCenter);


% THIS IS IN RUN FUNC NOW - so we can manipulate it per trial - should we
% keep anything in dv.pa still???
% dv.pa.theta = 135; % angle from center of one object (in degrees) converted to radians in func (keep it this way?)
% r = sqrt((dv.pa.Dx)^2 + (dv.pa.Dy)^2); % degrees of visual angle (which r will be the x and y dimensions)
% dv.pa.object1loc = [dv.pa.xCenter + r * cos(dv.pa.theta) dv.pa.yCenter + r * sin(dv.pa.theta)];
% dv.pa.object2loc = [dv.pa.xCenter + r * cos(dv.pa.theta+pi) dv.pa.yCenter + r * sin(dv.pa.theta+pi)];

%% Timing

% Task Events - (was in ms, keep it that way?)
dv.pa.sceneTime = 1;
dv.pa.objectSceneTime = 2;
dv.pa.probeCueTime = .5;
dv.pa.delayTime = 5;
dv.pa.probeTime = 2;

% Reward time: is time that solenoid is opened for. set to 100 miliseconds  
dv.pa.rewardTime = .1;
dv.pa.rewardWait = 0; 
dv.pa.breakFixPenalty = 2;
dv.pa.jitterSize = .5;

% fixation
dv.pa.preTrial = .5;       
dv.pa.fixWait = 4;        
dv.pa.fixHold = 1;
dv.pa.fixationXY = [0 0]; % degrees

% window sizes 
dv.pa.fpWin = [1.5 1.5]; % x,y radius
dv.pa.winScale = .2; % scale targwin with eccentricity

% dot sizes for drawing
dv.pa.eyeW      = 8;    % eye indicator width in pixels
dv.pa.fixdotW   = 8;    % width of the fixation dot
dv.pa.targdotW  = 8;    % width of the target dot
dv.pa.cursorW   = 8;   % cursor width in pixels



%% States - dv.states
dv.states.START     = 1;
dv.states.FPON      = 2;
dv.states.FPHOLD    = 3;

dv.states.VIRGINTRIAL = 4;

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


