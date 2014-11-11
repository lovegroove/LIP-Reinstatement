function [PDS,dv] = runLIPreinTrial(PDS,dv)
% This function is a single trial, trials are interated through the wrapper function runPLDAPS
% PDS is a struct containing your data
% dv is a struct containing your parameters (initialized in another function called lipreinstatement) and does not change
% dv.trial can be updated per trial to change paramters on each trial


dv = defaultTrialVariables(dv); % setup default trial struct (need this?)


%% Preallocation 
%-------------------------------------------------------------------------%
PDS.data.pairs = cell(dv.finish,4); % For saving correct objects and their placement
PDS.data.objectLocs = cell(dv.finish,1);
PDS.breakRestart = cell(dv.finish,1);
PDS.nBreaks = cell(dv.finish,1);
dv.trial.objectLocs = {}; % not sure of dimensions right now

% preallocate data aquisition variables
flipTimes       = zeros(1e4,2);
loopTimes       = zeros(1e4,2);
photodiodeTimes = zeros(1e4,2);
eyepos          = zeros(1e4,4);     % preallocate eye position data

% ISSUE - where does inputSamplingRate come from? when using eyetracking
% funcs? perhaps add a boolean for use only when using eyetracking
if dv.pass == 0
dv.disp.inputSamplingRate = 240;
eyeTimeStep     = 1/dv.disp.inputSamplingRate;
end

% add jitter
dv.trial.preTrial   = dv.pa.preTrial + rand*dv.pa.jitterSize;
dv.trial.fpWait     = dv.pa.fixWait;
dv.trial.fpHold     = dv.pa.fixHold + rand*dv.pa.jitterSize;
dv.trial.breakFixPenalty = dv.pa.breakFixPenalty;
dv.trial.graceTime = dv.pa.graceTime; % in case we want to add jitter

%-------------------------------------------------------------------------%
%Trial Parameters

dv.trial.breakRestart = 0; % (binary option) Restart trial when fixation is broken rather than go to next trail, so we don't miss pairings

dv.trial.nBreaks = 0; % Initialize

dv.trial.virgin = 1; % load images/make textures the first time through each trial
% or use dv.trial.state == dv.trial.states.VIRGINTRIAL******

% Initialize state times each trial (using different timers for these
% states on purpose right now)
dv.trial.showCueTime = 0;
dv.trial.showPairTime = 0;
dv.trial.delayTime = 0;
dv.trial.showProbeTime = 0;

%-------------------------------------------------------------------------%
% Stimulus (object) angle and distance from center (fixation pt)**********
dv.trial.theta = 135; % angle from center of one object (in degrees) converted to radians in func (can randomize this per trial if desired)
dv = stimLoc(dv);
% Make this a seperate stimulus lcoation support function???
% r = sqrt((dv.pa.Dx)^2 + (dv.pa.Dy)^2); % degrees of visual angle (which r will be the x and y dimensions)
% dv.trial.object1loc = [dv.pa.xCenter + r * cos(dv.trial.theta) dv.pa.yCenter + r * sin(dv.trial.theta)];
% dv.trial.object2loc = [dv.pa.xCenter + r * cos(dv.trial.theta+pi) dv.pa.yCenter + r * sin(dv.trial.theta+pi)];


%-------------------------------------------------------------------------%
%************************************ need this stuff????****************
%%% Degrees to pixels %%%
% We input positions on the screen in degrees of visual angle, but
% PsychToolbox takes calls to pixel locations. dv.disp.ppd is the pixels
% per degree of this monitor, given the viewing distance. It is important
% to set this up correctly in the setupPLADPS script and it requires
% Init_StereoDispPI to be implemented properly.
dv.trial.fpWin = dv.pa.fpWin*dv.disp.ppd;
% set fixation pixels
if ~isfield(dv.pa, 'fixationXY')
    dv.trial.fixXY = zeros(1,2); % fixation set to 0,0
else
    dv.trial.fixXY = dv.disp.ppd*dv.pa.fixationXY;
end

%% INITIALIZE STUFF
%-------------------------------------------------------------------------%
%%% setup PsychPortAudio %%%
%-------------------------------------------------------------------------%
pdsAudioClearBuffer(dv)
%-------------------------------------------------------------------------%
%%% Initalize Keyboard %%%
%-------------------------------------------------------------------------%
dv = pdsKeyboardClearBuffer(dv);
%-------------------------------------------------------------------------%
%%% Spike server %%%
%-------------------------------------------------------------------------%
[dv,spikes] = pdsSpikeServerGetSpikes(dv);
%-------------------------------------------------------------------------%
%%% Eyelink Toolbox Setup %%%
%-------------------------------------------------------------------------%
% preallocate for all eye samples and event data from the eyelink
dv = pdsEyelinkStartTrial(dv);
%-------------------------------------------------------------------------%
%%% Start Datapixx Analog Recording
%-------------------------------------------------------------------------%
if isfield(dv, 'dp')
    dv.dp = pdsDatapixxAdcStart(dv.dp);
end

%% Start & Sync Timing
%-------------------------------------------------------------------------%
% record start of trial in Datapixx, Mac & Plexon
% each device has a separate clock

% At the beginning of each trial, strobe a unique number to the plexon
% through the Datapixx to identify each trial. Often the Stimulus display
% will be running for many trials before the recording begins so this lets
% the plexon rig sync up its first trial with whatever trial number is on
% for stimulus display.
% SYNC clocks
clocktime = fix(clock);
for ii = 1:6
    pdsDatapixxStrobe(clocktime(ii));
end
PDS.unique_number(dv.j,:) = clocktime;    % trial identifier
Screen('Flip', dv.disp.ptr, 0);           % re-sync flip
dv.trial.trstart = GetSecs;
PDS.timing.datapixxStartTime(dv.j) = Datapixx('Gettime');
pdsDatapixxFlipBit(dv.events.TRIALSTART);  % start of trial (Plexon)
% Ping the datapixx to get the datapixx clock time. The datapixx clock time
% is always running and really functions as the master clock for all
% experiments. This is used for allignment later.
% if dv.useDatapixxbool                   % start of trial (Datapixx)
%     PDS.datapixxstarttime(dv.j) = Datapixx('Gettime');
% end

if dv.pass == 0
PDS.timing.eyelinkStartTime(dv.j) = Eyelink('TrackerTime');
Eyelink('message', 'TRIALSTART');
end

HideCursor;
dv.trial.ttime  = GetSecs - dv.trial.trstart;
PDS.timing.syncTimeDuration(dv.j) = dv.trial.ttime;


% % Query the frame duration - MY EDITION (needed? what about datapixx?)
% *****************
ifi = Screen('GetFlipInterval', dv.disp.ptr);
vbl = Screen('Flip', dv.disp.ptr); %Initially synchronize with retrace, take base time in vbl

%% Main While Loop for the Trial******************************************
%---------------------------------------------------------------------%
while ~dv.trial.flagNextTrial && dv.quit == 0
% trial time
    dv.trial.ttime  = GetSecs - dv.trial.trstart;
    
    %%% get inputs %%%
    %---------------------------------------------------------------------%
    % get mouse/eyetracker data
    
    [dv.trial.cursorX,dv.trial.cursorY,dv.trial.isMouseButtonDown] = GetMouse;
    
    dv = pdsGetEyePosition(dv); % Eye position or Mouse Position if UseMouse is flagged (within func)
    
    %dv = pdsGetEyePosition(dv, updateQueue); % Need this update queue????????????????? 
    
    %%% TRIAL STATES %%% dv.trial.state
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
        
    % 0) VIRGIN TRIAL - load images/make textures for current trial
    dv = virginTrial(dv);
    
    % 1) TURN ON FIXATION DOT
    dv = turnOnFixationPoint(dv);
    
    % 2) WAIT FOR FIXATION
    dv = waitForFixation(dv);
    
    % 3) HOLD FIXATION
    dv = holdFixation(dv);
    
    % 4) SHOW CUE
    dv = showCue(dv);
    
    % 5) SHOW PAIR (Study Trial)
    dv = showPair(dv);
    
    % 5) DELAY (Test Trial)
    dv = delay(dv);
    
    % 6) SHOW PROBE (Test Trial)
    dv = showProbe(dv);
    
    % 7) Trial Outcome
    dv = trialComplete(dv);
    dv = breakFixation(dv);
    
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    
    if dv.pass == 0
        
        % Update Measurements - need this?????
        %---------------------------------------------------------------------%
        dv.trial.ttime = GetSecs - dv.trial.trstart;
        % update eye position (fast)
        if dv.trial.ttime > dv.trial.timeLastEyePos + eyeTimeStep
            eyepos(dv.trial.iSample,:) = [dv.trial.eyeX, dv.trial.eyeY, dv.trial.ttime, dv.trial.state];
            dv.trial.timeLastEyePos = dv.trial.ttime;
            dv.trial.iSample = dv.trial.iSample + 1;
        end
        
    end
    
    
    %% KEYBOARD
    %---------------------------------------------------------------------%
[dv.trial.pressedQ dv.trial.firstPressQ]=KbQueueCheck(); % fast
    
    if dv.trial.firstPressQ(dv.kb.rKey)       % R = cycled targets
        dv.trial.targUser = 0;
    elseif dv.trial.firstPressQ(dv.kb.mKey)
        pdsDatapixxAnalogOut(dv.pa.rewardTime)
        pdsDatapixxFlipBit(dv.events.REWARD);
        dv.trial.ttime = GetSecs - dv.trial.trstart;
        dv.trial.timeReward(dv.trial.iReward) = dv.trial.ttime;
        dv.trial.iReward = dv.trial.iReward + 1;
        PsychPortAudio('Start', dv.pa.sound.reward);
    elseif dv.trial.firstPressQ(dv.kb.uKey)   % U = user selected targets
        dv.trial.targUser = 1;
    elseif dv.trial.firstPressQ(dv.kb.pKey)   % P = pause
        dv.quit = 1;
        ShowCursor
        Screen('Flip', dv.disp.ptr);
    elseif dv.trial.firstPressQ(dv.kb.qKey) % Q = quit
        Screen('Flip', dv.disp.ptr);
        dv = pdsEyelinkFinish(dv);
        PDS.timing.timestamplog = PsychDataPixx('GetTimestampLog', 1);
        dv.quit = 2;
        ShowCursor
    end
    

    %% DRAWING
 %---------------------------------------------------------------------%
    
 % draw textures in respective state functions or here?
 
 % FLIP ( or should we be doing their weird datapixx stuff???)
 vbl = Screen('Flip', dv.disp.ptr, vbl + 0.5*ifi);
 
 %[ig ig flipTimes(dv.trial.iFrame,1) flipTimes(dv.trial.iFrame,2)] = Screen('Flip', dv.disp.ptr); % use this?  

 %-------------------------------------------------------------------------%
    loopTimes(dv.trial.iLoop, :) = [dv.trial.ttime dv.trial.state];
    dv.trial.iLoop = dv.trial.iLoop + 1;
    
end % END WHILE LOOP

Screen('Close'); % was necessary in prototype version because of drawing all the textures

PDS.timing.datapixxstoptime(dv.j) = Datapixx('GetTime');
PDS.timing.trialend(dv.j) = GetSecs- dv.trial.trstart;

if isfield(dv, 'dp')
    dv.dp = pdsDatapixxAdcStop(dv.dp);
end

ShowCursor;

%%% Flush KbQueue %%%
KbQueueStop();
KbQueueFlush();

dv.goodtrial = dv.goodtrial + dv.trial.good;

%% Build PDS struct (your data and relevant timings and such) for saving %%
%-------------------------------------------------------------------------%
PDS.trialnumber(dv.j)        = dv.j;
PDS.goodtrial(dv.j)          = dv.trial.good;
PDS.fpXY(dv.j,:)             = (1/dv.disp.ppd)*dv.trial.fixXY;

PDS.timing.fpon(dv.j,:)      = [dv.trial.timeFpOn dv.trial.frameFpOn];
PDS.timing.fpentered(dv.j)   = dv.trial.timeFpEntered;
PDS.timing.fpoff(dv.j,:)     = [dv.trial.timeFpOff     dv.trial.frameFpOff];

PDS.timing.reward{dv.j}      = dv.trial.timeReward(~isnan(dv.trial.timeReward));

PDS.timing.breakfix(dv.j)    = dv.trial.timeBreakFix;

% system timing
PDS.timing.photodiodeFlips{dv.j} = photodiodeTimes(1:dv.trial.iPhotodiode -1,:);
PDS.timing.pldapsLoopTimes{dv.j}    = loopTimes(1:dv.trial.iLoop-1,:);
PDS.timing.ptbTrialStart(dv.j)      = dv.trial.trstart;
PDS.timing.ptbFliptimes{dv.j}       = flipTimes(flipTimes(:,1)~=0,:);
% PDS.timing.stimRefreshTime{dv.j} =
% stimRefreshTime(~isnan(stimRefreshTime)); % need this?

if isfield(dv, 'dp')
    analogTime = linspace(dv.dp.adctstart, dv.dp.adctend, size(dv.dp.bufferData,2));
    PDS.data.datapixxAnalog{dv.j} = [analogTime(:) dv.dp.bufferData'];
end

if dv.pass == 0
    %Eyelink
    PDS.data.eyeposLoop{dv.j} = eyepos;
    PDS.data.eyelinkSampleBuffer{dv.j}  =  dv.el.sampleBuffer(:,~isnan(dv.el.sampleBuffer(1,:)));
    PDS.data.eyelinkEventBuffer{dv.j}   =  dv.el.eventBuffer(:,~isnan(dv.el.eventBuffer(1,:)));
    PDS.data.eyelinkQueueStruct = eye;
    
    PDS.eyepos{dv.j}         = eyepos(1:dv.trial.iSample+1,:);  % remove extra    
end

fprintf(' %.0f/%.0f, %.2f good.\r', sum(PDS.goodtrial), length(PDS.goodtrial), (sum(PDS.goodtrial)/length(PDS.goodtrial)))

% SAVE pairings per trial and correct and foil object locations
PDS.data.pairs{dv.j} = dv.pairOrder(dv.j,:);   % is this stuff right?
PDS.data.objectLocs{dv.j} = dv.trial.objectLocs; 

PDS.breakRestart{dv.j} = dv.trial.breakRestart;
PDS.nBreaks{dv.j} = dv.trial.nBreaks;


%-------------------------------------------------------------------------%
%%% INLINE FUNCTIONS
%-------------------------------------------------------------------------%
    function held = fixationHeld(dv)
        %         held = squarewindow(dv.pass,dv.disp.ctr(1:2)+dv.trial.fixXY-[dv.trial.eyeX dv.trial.eyeY],dv.trial.winW,dv.trial.winH);
        %         held = squarewindow(dv.pass,dv.disp.ctr(1:2)+dv.trial.fixXY-[dv.trial.eyeX dv.trial.eyeY],dv.trial.winW,dv.trial.winH);
        held = circlewindow(dv.pass,[dv.trial.eyeX dv.trial.eyeY]-(dv.trial.fixXY+dv.disp.ctr(1:2)),dv.trial.fpWin(1),dv.trial.fpWin(2));
    end

%%% STATE FUNCTIONS %%%
%
% TURN ON FIXATION POINT
%---------------------------------------------------------------------%
    function dv = turnOnFixationPoint(dv)
        % BEFORE FIXATION DOT
        if  dv.trial.state == dv.states.START && dv.trial.ttime  > dv.trial.preTrial
            dv.trial.colorFixDot  = dv.disp.clut.targetgood;
            dv.trial.colorFixWindow = dv.disp.clut.bg;
            
            %%% TURN ON FIXATION DOT %%%
            pdsDatapixxFlipBit(dv.events.FIXATION) % fp1 ON
            dv.trial.timeFpOn = GetSecs - dv.trial.trstart;   % time state was entered
            dv.trial.frameFpOn = dv.trial.iFrame; % keep track of the frame that this will be displayed on
            dv.trial.state = dv.states.FPON;
            
            if dv.sound
                % Start of the trial with a cue sound.
                PsychPortAudio('Start', dv.pa.sound.cue);
            end
        end
    end
% WAIT FOR FIXATION
%---------------------------------------------------------------------%
    function dv = waitForFixation(dv)
        dv.trial.fpWin = dv.pa.fpWin*dv.disp.ppd/2;
        % WAITING FOR SUBJECT FIXATION (fp1)
        if  dv.trial.state == dv.states.FPON
            if dv.trial.ttime  < (dv.trial.preTrial+dv.trial.fpWait) && fixationHeld(dv)
                dv.trial.colorFixDot = dv.disp.clut.targetnull;
                dv.trial.colorFixWindow = dv.disp.clut.window;
                dv.trial.timeFpEntered = GetSecs - dv.trial.trstart;
                pdsDatapixxFlipBit(dv.events.FIXATION)
                dv.trial.state = dv.states.FPHOLD;
            elseif dv.trial.ttime  > (dv.trial.preTrial+dv.trial.fpWait) % Was inflexible. TBC: dv.pa.fixwaitstop
                pdsDatapixxFlipBit(dv.events.BREAKFIX)
                dv.trial.timeBreakFix = GetSecs - dv.trial.trstart;
                dv.trial.state = dv.states.BREAKFIX;
            end
        end
        
        
    end

% HOLD FIXATION FOR APPROPRIATE DURATION
%---------------------------------------------------------------------%
    function dv = holdFixation(dv)
        dv.trial.fpWin = dv.pa.fpWin*dv.disp.ppd;
        % check if fixation is held
        if dv.trial.state == dv.states.FPHOLD
            if dv.trial.ttime > dv.trial.timeFpEntered && fixationHeld(dv) % was fixPtOffset + timeFpEntered but I don't think that applies to mine
                dv.trial.colorFixDot    = dv.disp.clut.bg;
                dv.trial.colorFixWindow = dv.disp.clut.bg;
                pdsDatapixxFlipBit(dv.events.FIXATION) % fixation cross
                dv.trial.ttime      = GetSecs - dv.trial.trstart;
                dv.trial.timeFpOff  = dv.trial.ttime;
                dv.trial.state      = dv.states.SHOWCUE;  % switching state to show cue (scene)
            elseif dv.trial.ttime < dv.trial.timeFpEntered && ~fixationHeld(dv)
                dv.trial.colorFixDot    = dv.disp.clut.bg;
                dv.trial.colorFixWindow = dv.disp.clut.bg;
                pdsDatapixxFlipBit(dv.events.BREAKFIX)
                dv.trial.timeBreakFix = GetSecs - dv.trial.trstart;
                dv.trial.state = dv.states.BREAKFIX;
            end
        end
    end
    
    % TRIAL COMPLETE -- GIVE REWARD  -
    % ********************************************** stuff in here we don't
    % need/want???
%---------------------------------------------------------------------%
    function dv = trialComplete(dv)
        if dv.trial.state == dv.states.TRIALCOMPLETE
%             dv.trial.colorTarget1Dot    = dv.disp.clut.bg;            % target color
%             dv.trial.colorTarget2Dot    = dv.disp.clut.bg;            % target color
            dv.trial.colorFixWindow     = dv.disp.clut.bg;           % fixation window color
%             dv.trial.colorTarget1Window = dv.disp.clut.bg;
%             dv.trial.colorTarget2Window = dv.disp.clut.bg;
            dv.trial.good = 1;
            if dv.pa.freeOn == 0
                if dv.trial.targ1Chosen == dv.trial.targ1Correct
                    dv.trial.correct = 1;
                    soundHandle = dv.pa.sound.reward;
                else
                    dv.trial.correct = 0;
                    soundHandle = dv.pa.sound.incorrect;
                end
                % if free trial, give reward independent of choice:
            elseif dv.pa.freeOn == 1
                dv.trial.correct = 1;
                soundHandle = dv.pa.sound.reward;
            end
            
            if dv.trial.ttime > dv.trial.timeComplete + .2
                
                if dv.trial.flagBuzzer
                    pdsDatapixxFlipBit(dv.events.REWARD);
                    PsychPortAudio('Start', soundHandle, 1, [], [], GetSecs + .1);
                    dv.trial.ttime = GetSecs - dv.trial.trstart;
                    dv.trial.timeReward(dv.trial.iReward) = dv.trial.ttime;
                    dv.trial.iReward = dv.trial.iReward + 1;
                    if dv.trial.correct
                        pdsDatapixxAnalogOut(dv.pa.rewardTime)
                    end
                    dv.trial.flagBuzzer = 0;
                end
            end
            
            if dv.trial.ttime > dv.trial.timeComplete + 1
                pdsDatapixxFlipBit(dv.events.TRIALEND);
                dv.trial.flagNextTrial = true;
            end
            
        end
        
        
        
    end

% BREAK FIXATION
%---------------------------------------------------------------------%
    function dv = breakFixation(dv)
        if dv.trial.state == dv.states.BREAKFIX
            % turn off stimulus
            dv.trial.colorFixDot        = dv.disp.clut.bg;            % fixation point 1 color
            dv.trial.colorTarget1Dot    = dv.disp.clut.bg;            % target color
            dv.trial.colorTarget2Dot    = dv.disp.clut.bg;            % target color
            dv.trial.colorFixWindow     = dv.disp.clut.bg;           % fixation window color
            dv.trial.colorTarget1Window = dv.disp.clut.bg;
            dv.trial.colorTarget2Window = dv.disp.clut.bg;
            dv.trial.Gpars(4,:) = 0; % set gabors to on contrast
            dv.trial.Mpars(4,:) = 0; % set gabors to on contrast
            dv.trial.good = 0;
            dv.trial.flagMotionOn = 2;
            dv.trial.targOn = 2;
            audiohandle = dv.pa.sound.breakfix;
            if dv.trial.flagBuzzer
                PsychPortAudio('Start', audiohandle, 1, [], [], GetSecs + .1);
                dv.trial.flagBuzzer = 0;
            end
            
            if dv.trial.ttime > dv.trial.timeBreakFix + dv.trial.breakFixPenalty
                if dv.trial.breakRestart
                    dv.trial.state = dv.states.START; % option to restart trial to fix below, and probably cause other issues***
                    dv.trial.nBreaks = dv.trial.nBreaks + 1;
                else
                dv.trial.flagNextTrial = true; % currently breaking fixation increments the trial so you would miss the presentation of certain pairs in my paradigm
                end
                
            end
        end
    end


% VIRGIN TRIAL - Load Images and Make Textures for current trial
%---------------------------------------------------------------------%
    function dv = virginTrial(dv)
        if dv.trial.virgin % change to virgin trial state?????? for consistancy but this might make more trouble
            %         if dv.trial.state == dv.trial.states.VIRGINTRIAL
            %         end
            
            % Load in Image files
            dv.trial.objectImageFile = [dv.fileInfo.objectPath dv.pairOrder{dv.j,1}];
            dv.trial.sceneImageFile = [dv.fileInfo.scenePath dv.pairOrder{dv.j,2}];
            
            % Background
            dv.trial.sceneImage = imread(dv.trial.sceneImageFile, dv.fileInfo.sceneFormat); % Read in image
            dv.trial.sceneImage = imresize(dv.trial.sceneImage, [dv.disp.winRect(4) dv.disp.winRect(3)]); % Resize Image to screen - [y,x] weird
            dv.trial.sceneImageTexture = Screen('MakeTexture', dv.disp.ptr, dv.trial.sceneImage); % Create Texture
            %---------------------------------------------------------------------%
            % For STUDY Trial Type
            if strcmp(dv.trialType, 'study')
                
                % Object - make texture and assign location
                dv.trial.objectImage = imread(dv.trial.objectImageFile, dv.fileInfo.objectFormat); % Read in image
                dv.trial.objectImageTexture = Screen('MakeTexture', dv.disp.ptr, dv.trial.objectImage); % Create Texture
                [s1, s2, s3] = size(dv.trial.objectImage);
                dv.trial.baseRect = [0 0 s1 s2];
                dv.trial.destRect = CenterRectOnPointd(dv.trial.baseRect .* dv.pa.objectSize, dv.pa.xCenter, dv.pa.yCenter); %Set size and location (find a better way to do size manipulations)
                
            %---------------------------------------------------------------------%
                % For TEST Trial Type
            elseif strcmp(dv.trialType, 'test')
                
                % Randomly select foil and place object images randomly - can choose number of objects here
                dv.trial.foilObject = dv.pairOrder{randsample(length(dv.pairOrder),1) ,1};
                
                % Also reroll here just in case you happen to pick the same object
                while strcmp(dv.trial.foilObject, dv.pairOrder{dv.j,1})
                    dv.trial.foilObject = dv.pairOrder{randsample(length(dv.pairOrder),1) ,1};
                end
                % Place objects randomly
                dv.trial.foilObjectImageFile = [dv.fileInfo.objectPath dv.trial.foilObject];
                dv.trial.placedObjects = cell(2,1);
                dv.trial.placedObjects{1} = dv.trial.objectImageFile;
                dv.trial.placedObjects{2} = dv.trial.foilObjectImageFile;
                [dv.trial.object1, dv.trial.object2] = dv.trial.placedObjects{randsample(2,2)};
                
                % Object 1 - make texture and assign location
                dv.trial.object1Image = imread(dv.trial.object1, dv.fileInfo.objectFormat); % Read in image
                dv.trial.object1ImageTexture = Screen('MakeTexture', dv.disp.ptr, dv.trial.object1Image); % Create Texture
                [s1, s2, s3] = size(dv.trial.object1Image);
                dv.trial.baseRect = [0 0 s1 s2];
                dv.trial.destRect1 = CenterRectOnPointd(dv.trial.baseRect .* dv.pa.objectSize, dv.trial.object1loc(1), dv.trial.object1loc(2)); %Set size and location
                
                % Object 2 - make texture and assign location
                dv.trial.object2Image = imread(dv.trial.object2, dv.fileInfo.objectFormat); % Read in image
                dv.trial.object2ImageTexture = Screen('MakeTexture', dv.disp.ptr, dv.trial.object2Image); % Create Texture
                [s1, s2, s3] = size(dv.trial.object2Image);
                dv.trial.baseRect = [0 0 s1 s2];
                dv.trial.destRect2 = CenterRectOnPointd(dv.trial.baseRect .* dv.pa.objectSize, dv.trial.object2loc(1), dv.trial.object2loc(2)); %Set size and location
                
                % Save correct object and location
                dv.trial.objectLocs = {dv.trial.placedObjects{1}, dv.trial.placedObjects{2}; dv.trial.destRect1, dv.trial.destRect2};
                
            end % end trial type if
            
            dv.trial.virgin = 0;
        end % End virgin if 
    end % END virgin func

% SHOW CUE
%---------------------------------------------------------------------%
    function dv = showCue(dv)
        % shows Scene Cue in both study and test trials
        
        if dv.trial.state == dv.states.SHOWCUE 
            
            % using internal timers for each state to facilitate modularity
            % and future changes (might not be the best choice
            if dv.trial.showCueTime == 0 % reset timer first time through
                dv.trial.showCueTime = 0;
            end
            dv.trial.showCueTime = GetSecs - dv.trial.showCueTime;
            
            %%%%%% Scene Alone (1s default)
            Screen('DrawTexture',dv.disp.ptr,dv.trial.sceneImageTexture); % Draw Texture (same texture variable for both trial types?)
            
            if dv.trial.showCueTime >= dv.pa.sceneTime
                
                if strcmp(dv.trialType, 'study')
                    dv.trial.state = dv.states.SHOWPAIR;
                elseif strcmp(dv.trialType, 'test')
                    dv.trial.state = dv.states.DELAY;
                end
                
            end
            
        end %end if trial state
    end % end func

% SHOW PAIR
%---------------------------------------------------------------------%

    function dv = showPair(dv)
        if dv.trial.state == dv.states.SHOWPAIR 
            % Just for study trials
            
            if dv.trial.showPairTime == 0 % reset timer first time through
                dv.trial.showPairTime = 0; 
            end
            dv.trial.showPairTime = GetSecs - dv.trial.showPairTime;
            
            %%%%%% Scene and Object (2s default)
            % Scene
            Screen('DrawTexture',dv.disp.ptr,dv.trial.sceneImageTexture);
            % Object
            Screen('DrawTexture',dv.disp.ptr,dv.trial.objectImageTexture,[],dv.trial.destRect); % Draw Texture
            
            if dv.trial.showPairTime >= dv.pa.objectSceneTime
                dv.trial.state = dv.states.TRIALCOMPLETE;
            end
            
        end %end if trial state
    end % end func

 % DELAY - this is going to be a pain in the ass to deal with breaking
 % fixation etc.!!!!!!!!!!!!!!!!!! and what about brief fix pt cue after?
 %---------------------------------------------------------------------%
    function dv = delay(dv)
        if dv.trial.state == dv.states.DELAY
            if dv.trial.delayTime == 0 % reset timer
                dv.trial.delayTime = 0; 
            end
            dv.trial.delayTime = GetSecs - dv.trial.delayTime;
            
            if fixationHeld(dv) || dv.trial.delayTime < dv.trial.graceTime % need a grace time for them to look back at the center! (does this work alright?)
                
                if dv.trial.delayTime <= dv.pa.delayTime
                    %%%% Delay (5s default)
                    Screen('FillRect', dv.disp.ptr, dv.pa.centeredDelayRect); % what to do about color here?
                    
                elseif dv.trial.delayTime > dv.pa.delayTime % briefly represent fixation pt to cue onset of memory probe
                    
                    %%% TURN ON FIXATION DOT %%% - does this work??????????? is this because of weird drawing all the time as bg color
                    pdsDatapixxFlipBit(dv.events.FIXATION) % fp1 ON
                    
                elseif dv.trial.delayTime >= dv.pa.delayTime + dv.pa.probeCueTime
                    dv.trial.state = dv.states.SHOWPROBE;
                end
                
            else
                dv.trial.state = dv.states.BREAKFIX;
            end
        end
    end
    
% SHOW PROBE
%---------------------------------------------------------------------%
    function dv = showProbe(dv)
        if dv.trial.state == dv.states.SHOWPROBE
            % just for test trials
            
            if dv.trial.showProbeTime == 0 % reset timer
                dv.trial.showProbeTime = 0; 
            end
            
            dv.trial.showProbeTime = GetSecs - dv.trial.showProbeTime;
        
                %%%%%% Scene and 2 Objects (2s)
                % ******** Scene ********
                Screen('DrawTexture',dv.disp.ptr,sceneImageTexture);
                
                % ******** Object 1 ********
                Screen('DrawTexture',dv.disp.ptr,dv.trial.object1ImageTexture,[],dv.trial.destRect1); % Draw Texture
                
                % ******** Object 2 ********
                Screen('DrawTexture',dv.disp.ptr,dv.trial.object2ImageTexture,[],dv.trial.destRect2); % Draw Texture
        end
        
        if dv.trial.showProbeTime >= dv.pa.probeTime
            dv.trial.state = dv.states.TRIALCOMPLETE;
        end
        
    end





end % end master function