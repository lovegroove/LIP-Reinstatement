function [PDS,dv] = runLIPreinTrial(PDS,dv)
% This function is a single trial, trials are interated through the wrapper
% function runPLDAPScv (cv = color version)
% PDS is a struct containing your data
% dv is a struct containing your parameters (initialized in another function called liprein) and does not change
% dv.trial can be updated per trial to change paramters on each trial


dv = defaultTrialVariables(dv); % setup default trial struct


%% Preallocation 
%-------------------------------------------------------------------------%
dv.trial.objectLocs = {}; % not sure of dimensions right now

% preallocate data aquisition variables
flipTimes       = zeros(1e4,2);
loopTimes       = zeros(1e4,2);
photodiodeTimes = zeros(1e4,2);
eyepos          = zeros(1e4,4);     % preallocate eye position data

if dv.pass == 0
dv.disp.inputSamplingRate = 240;
eyeTimeStep     = 1/dv.disp.inputSamplingRate;
end

% add jitter
dv.trial.preTrial   = dv.pa.preTrial + rand*dv.pa.jitterSize;
dv.trial.fixPtOffset = dv.pa.fixPtOffset; %+ rand*dv.pa.jitterSize;
dv.trial.fpWait     = dv.pa.fixWait;
dv.trial.fpHold     = dv.pa.fixHold + rand*dv.pa.jitterSize;
dv.trial.breakFixPenalty = dv.pa.breakFixPenalty;
dv.trial.graceTime = dv.pa.graceTime; % in case we want to add jitter

%-------------------------------------------------------------------------%
%Trial Parameters

% Special Break fixation Options so we don't "miss" pairs being presented
dv.trial.breakRestart = 0; % (binary option) Restart trial when fixation is broken rather than go to next trail, so we don't miss pairings
dv.trial.nBreaks = 0; % Initialize

% Set flags
dv.trial.fixFlagOn = 0;
dv.trial.fixDotOn = 0;
dv.trial.delayRectOn = 0;
dv.trial.virgin = 1; % boolean, load images/make textures the first time through each trial
dv.trial.showCueFlag = 1;
dv.trial.showPairFlag = 1;
dv.trial.delayFlag = 1;
dv.trial.showProbeFlag = 1;
dv.trial.waitFPDelayFlag = 1;
dv.trial.holdFixDelayFlag = 1;
dv.trial.eyePosi = 0; %initialize counter
dv.trial.switchTrialType = 0;
%-------------------------------------------------------------------------%
% Stimulus (object) angle and distance from center (fixation pt), foil
% object placed 180 deg. from correct object

%randomize each trial if desired
dv.trial.theta = dv.pa.stimThetas(randi(length(dv.pa.stimThetas)));

%dv.trial.theta = 100; % angle from center of one object (in degrees) converted to radians in func
dv = stimLoc(dv);

% Delay timing if random - Overwrites dv.pa.delayTime 
if dv.useRandomDelay
dv.pa.delayTime = randi(5)+.25-1; % current range 0.25-4.25 because it will dv.pa.probeCueTime will add .75 
end
%-------------------------------------------------------------------------%
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


% % Query the frame duration - needed?
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
    
    % Alternative Strict Delay (Test Trial)
    dv = waitForFixationDelay(dv);
    dv = holdFixationDelay(dv);
    
    % 6) SHOW PROBE (Test Trial)
    dv = showProbe(dv);
    
    % 7) Trial Outcome
    dv = trialComplete(dv);
    dv = breakFixation(dv);
    
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    
    if dv.pass == 0
        
        % Update Measurements
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
    elseif dv.trial.firstPressQ(dv.kb.tKey) % T = switch trial type
        dv.quit = 1;
        ShowCursor
        Screen('Flip', dv.disp.ptr);
        dv.trial.switchTrialType = 1;
        disp('Trial type switched. Type return to start next trial...')
    end
    

    %% DRAWING
 %---------------------------------------------------------------------%
 % Show fixation cross if desired 
 if dv.trial.fixFlagOn
     Screen('DrawLines', dv.disp.ptr, dv.pa.allCoords, dv.pa.lineWidthPix, dv.pa.fixCrossColor,  [dv.pa.xCenter, dv.pa.yCenter], 2)
 end
 
 % Show red dot on fixation cross
 if dv.trial.fixDotOn
     Screen('DrawDots', dv.disp.ptr, [dv.pa.xCenter, dv.pa.yCenter], 10, [1,0,0], [], 2)
 end
 
 % Show Mouse pointer if desired 
 if dv.showMouse
     Screen('DrawDots', dv.disp.ptr, [dv.trial.cursorX,dv.trial.cursorY], 10, [1,1,1], [], 2);
 end
    
 % FLIP
 vbl = Screen('Flip', dv.disp.ptr, vbl + 0.5*ifi);
 
 %-------------------------------------------------------------------------%
    loopTimes(dv.trial.iLoop, :) = [dv.trial.ttime dv.trial.state];
    dv.trial.iLoop = dv.trial.iLoop + 1;
    
end % END WHILE LOOP

Screen('Close'); % was necessary in prototype version because of drawing all the textures, need to close offscreen windows

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


% Saving timing states of interest and trial type for each trial
if PDS.goodtrial(dv.j) == 1
    PDS.timing.timeShowCueStart{dv.j} = dv.trial.timeShowCueStart;
    if strcmp(dv.trialType, 'study')
        PDS.trialType{dv.j} = dv.trialType;
        PDS.timing.timeShowPairStart{dv.j} = dv.trial.timeShowPairStart;
    elseif strcmp(dv.trialType, 'test')
        PDS.trialType{dv.j} = dv.trialType;
        PDS.timing.timeDelayStart{dv.j} = dv.trial.timeDelayStart;
        PDS.timing.timeShowProbeStart{dv.j} = dv.trial.timeShowProbeStart;
        PDS.timing.randomDelayLength{dv.j} = dv.pa.delayTime + dv.pa.probeCueTime;
    end
else
    PDS.timing.timeShowCueStart{dv.j} = 0;
    if strcmp(dv.trialType, 'study')
        PDS.trialType{dv.j} = dv.trialType;
        PDS.timing.timeShowPairStart{dv.j} = 0;
    elseif strcmp(dv.trialType, 'test')
        PDS.trialType{dv.j} = dv.trialType;
        PDS.timing.timeDelayStart{dv.j} = 0;
        PDS.timing.timeShowProbeStart{dv.j} = 0;
    end
end

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

% Save quick eye measures
if strcmp(dv.trialType, 'test')
PDS.data.eyeLocProbe{dv.j} = dv.trial.eyeLoc;
PDS.data.correctObject(dv.j) = dv.trial.correctObject;
end

%fprintf(' %.0f/%.0f, %.2f good.\r', sum(PDS.goodtrial), length(PDS.goodtrial), (sum(PDS.goodtrial)/length(PDS.goodtrial)))

% SAVE pairings per trial and correct and foil object locations
PDS.data.pairs{dv.j} = dv.pairOrder(dv.j,:);  
PDS.data.objectLocs{dv.j} = dv.trial.objectLocs; 

% depreciated
% PDS.timing.breakRestart{dv.j} = dv.trial.breakRestart;
% PDS.timing.nBreaks{dv.j} = dv.trial.nBreaks;

% For a Single Session
if dv.singleSession && dv.j == dv.pa.singleSessionStudy && strcmp(dv.trialType,'study') %current number of pairs in a block, make this a variable later
    disp('Study Session finished. Test time!')
    dv.quit = 1;
    ShowCursor
    
elseif dv.singleSession && dv.j == dv.pa.singleSessionTest && strcmp(dv.trialType,'test')
    disp('Test Session finished. Thank you!')
    dv.quit = 2;
    ShowCursor
end

% For naive human subject
if dv.naiveSubj && dv.j == dv.pa.singleSessionStudy
    disp('Session 1 finished. Please see experimenter to start session 2.')
    dv.trial.switchTrialType = 1;
    dv.quit = 1; % could leave this out and make it continuous 
    ShowCursor
    elseif dv.naiveSubj && dv.j == dv.pa.singleSessionTest 
    disp('Session 2 finished. Thank you!')
    dv.quit = 2;
    ShowCursor
end

% Quick trial type switching when t is pressed, it automatically toggles it
if dv.trial.switchTrialType
    if strcmp(dv.trialType, 'study')
        dv.trialType = 'test';
    elseif strcmp(dv.trialType, 'test')
        dv.trialType = 'study';
    end
    dv.trial.switchTrialType = 0; % probably not needed since set at the beginning of the trial as well
end

%-------------------------------------------------------------------------%
%%% INLINE FUNCTIONS
%-------------------------------------------------------------------------%
    function held = fixationHeld(dv)
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
%             dv.trial.colorFixDot  = dv.disp.clut.red;
%             dv.trial.colorFixWindow = dv.disp.clut.bg;
            
            dv.trial.fixFlagOn = 1;
            dv.trial.fixDotOn = 1;

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
            
            dv.trial.fixFlagOn = 1;
            dv.trial.fixDotOn = 1;
            
            if dv.trial.ttime  < (dv.trial.preTrial+dv.trial.fpWait) && fixationHeld(dv)
%                 dv.trial.colorFixDot = dv.disp.clut.targetnull;
%                 dv.trial.colorFixWindow = dv.disp.clut.window;
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
% WAIT FOR FIXATION FOR DELAY!!!
%---------------------------------------------------------------------%
function dv = waitForFixationDelay(dv)
        dv.trial.fpWin = dv.pa.fpWin*dv.disp.ppd/2;
        % WAITING FOR SUBJECT FIXATION (fp1)
        if  dv.trial.state == dv.states.FPONDELAY
            
            if dv.trial.waitFPDelayFlag
            dv.trial.ttime = GetSecs - dv.trial.trstart;
            dv.trial.waitFPDelayStart = dv.trial.ttime;
            dv.trial.waitFPDelayFlag = 0;
            end
            
            dv.trial.fixFlagOn = 1;
            dv.trial.fixDotOn = 1;
            
            if dv.trial.ttime  < (dv.trial.waitFPDelayStart+dv.pa.graceTime) && fixationHeld(dv)
                dv.trial.timeFpEntered = GetSecs - dv.trial.trstart; % probably overwrites FP times from beginning of the trial (do we care right now?)
                pdsDatapixxFlipBit(dv.events.FIXATION)
                dv.trial.state = dv.states.FPHOLDDELAY;
            elseif dv.trial.ttime  > (dv.trial.waitFPDelayStart + dv.pa.graceTime)
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
            
            dv.trial.fixFlagOn = 1;
            dv.trial.fixDotOn = 1;
            
            if dv.trial.ttime > dv.trial.timeFpEntered + dv.trial.fixPtOffset && fixationHeld(dv)
%                 dv.trial.colorFixDot    = dv.disp.clut.bg;
%                 dv.trial.colorFixWindow = dv.disp.clut.bg;
                dv.trial.fixFlagOn = 0;
                dv.trial.fixDotOn = 0;
                
                pdsDatapixxFlipBit(dv.events.FIXATION) % fixation cross
                dv.trial.ttime      = GetSecs - dv.trial.trstart;
                dv.trial.timeFpOff  = dv.trial.ttime;
                dv.trial.state      = dv.states.SHOWCUE;  % switching state to show cue (scene)
            elseif dv.trial.ttime < dv.trial.timeFpEntered && ~fixationHeld(dv)
%                 dv.trial.colorFixDot    = dv.disp.clut.bg;
%                 dv.trial.colorFixWindow = dv.disp.clut.bg;
                dv.trial.fixFlagOn = 0;
                dv.trial.fixDotOn = 0;
                pdsDatapixxFlipBit(dv.events.BREAKFIX)
                dv.trial.timeBreakFix = GetSecs - dv.trial.trstart;
                dv.trial.state = dv.states.BREAKFIX;
            end
        end
    end
    
    
    % HOLD FIXATION FOR APPROPRIATE DURATION DURING DELAY!!!!
%---------------------------------------------------------------------%
    function dv = holdFixationDelay(dv)
        dv.trial.fpWin = dv.pa.fpWin*dv.disp.ppd;
        % check if fixation is held
        if dv.trial.state ==  dv.states.FPHOLDDELAY
            
            if dv.trial.holdFixDelayFlag
            dv.trial.ttime = GetSecs - dv.trial.trstart;
            dv.trial.holdFixDelayStart = dv.trial.ttime;
            dv.trial.holdFixDelayFlag = 0;
            end
            
            dv.trial.fixFlagOn = 1;
            dv.trial.fixDotOn = 1;
            
            if dv.trial.ttime > (dv.trial.holdFixDelayStart + dv.pa.delayTime + dv.pa.probeCueTime) && fixationHeld(dv)
                dv.trial.fixFlagOn = 0;
                dv.trial.fixDotOn = 0;
                pdsDatapixxFlipBit(dv.events.FIXATION) % fixation cross
                dv.trial.ttime      = GetSecs - dv.trial.trstart; 
                dv.trial.timeFpOff  = dv.trial.ttime; % Overwriting times on intitial fixation probably 
                dv.trial.state      = dv.states.SHOWPROBE;  
                
            elseif dv.trial.ttime < (dv.trial.holdFixDelayStart + dv.pa.delayTime + dv.pa.probeCueTime) && ~fixationHeld(dv)
                dv.trial.fixFlagOn = 0;
                dv.trial.fixDotOn = 0;
                pdsDatapixxFlipBit(dv.events.BREAKFIX)
                dv.trial.timeBreakFix = GetSecs - dv.trial.trstart;
                dv.trial.state = dv.states.BREAKFIX;
                
            elseif dv.trial.ttime > (dv.trial.holdFixDelayStart + dv.pa.delayTime) && fixationHeld(dv)
                 % Red dot off for visual cue
                dv.trial.fixDotOn = 0;
                dv.trial.fixFlagOn = 1;
            end
        end
    end
    
    
    % TRIAL COMPLETE -- GIVE REWARD  -
    % ********************************************** 
%---------------------------------------------------------------------%
    function dv = trialComplete(dv)
        if dv.trial.state == dv.states.TRIALCOMPLETE
            
            dv.trial.fixFlagOn = 0; %probably not needed
            dv.trial.fixDotOn = 0;
            
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
            
            dv.trial.fixFlagOn = 0;
            dv.trial.fixDotOn = 0;
            
            audiohandle = dv.pa.sound.breakfix;
            if dv.trial.flagBuzzer
                PsychPortAudio('Start', audiohandle, 1, [], [], GetSecs + .1);
                dv.trial.flagBuzzer = 0;
            end
            
            if dv.trial.ttime > dv.trial.timeBreakFix + dv.trial.breakFixPenalty
                if dv.trial.breakRestart
                    dv.trial.state = dv.states.START; % option to restart trial to fix below, and probably cause other issues***
                    dv.trial.nBreaks = dv.trial.nBreaks + 1;
                    dv.trial.trstart = GetSecs; % restart trial timer
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
                dv.trial.baseRect = [0 0 dv.disp.ppd .* dv.pa.objectSize dv.disp.ppd .* dv.pa.objectSize];
                dv.trial.destRect = CenterRectOnPointd(dv.trial.baseRect, dv.pa.xCenter, dv.pa.yCenter); %Set size and location (find a better way to do size manipulations)
                
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
                % Save correct
                if strcmp(dv.trial.object1,dv.trial.objectImageFile)
                    dv.trial.correctObject = 1;
                elseif strcmp(dv.trial.object2,dv.trial.objectImageFile)
                    dv.trial.correctObject = 2;
                end
                
                % Object 1 - make texture and assign location
                dv.trial.object1Image = imread(dv.trial.object1, dv.fileInfo.objectFormat); % Read in image
                dv.trial.object1ImageTexture = Screen('MakeTexture', dv.disp.ptr, dv.trial.object1Image); % Create Texture
                dv.trial.baseRect = [0 0 dv.disp.ppd .* dv.pa.objectSize dv.disp.ppd .* dv.pa.objectSize];
                dv.trial.destRect1 = CenterRectOnPointd(dv.trial.baseRect, dv.trial.object1loc(1), dv.trial.object1loc(2)); %Set size and location
                
                % Object 2 - make texture and assign location
                dv.trial.object2Image = imread(dv.trial.object2, dv.fileInfo.objectFormat); % Read in image
                dv.trial.object2ImageTexture = Screen('MakeTexture', dv.disp.ptr, dv.trial.object2Image); % Create Texture
                dv.trial.baseRect = [0 0 dv.disp.ppd .* dv.pa.objectSize dv.disp.ppd .* dv.pa.objectSize];
                dv.trial.destRect2 = CenterRectOnPointd(dv.trial.baseRect, dv.trial.object2loc(1), dv.trial.object2loc(2)); %Set size and location
                
                % Save correct object and location
                dv.trial.objectLocs = {dv.trial.placedObjects{1}, dv.trial.placedObjects{2}, dv.trial.destRect1, dv.trial.destRect2, dv.trial.correctObject}; % [correctObj, foilObj, loc1, loc2, correctLoc] locations are random by this point so we don't know which is which, added dv.trial.correctObject to a column here
                
            end % end trial type if
            
            dv.trial.virgin = 0;
        end % End virgin if 
    end % END virgin func

% SHOW CUE
%---------------------------------------------------------------------%
    function dv = showCue(dv)
        % shows Scene Cue in both study and test trials
        if dv.trial.state == dv.states.SHOWCUE
            
            dv.trial.fixFlagOn = 0;
            dv.trial.fixDotOn = 0;
            
            if dv.trial.showCueFlag
                dv.trial.ttime = GetSecs - dv.trial.trstart;
                dv.trial.timeShowCueStart = dv.trial.ttime;
                dv.trial.showCueFlag = 0;
            end
            
            %%%%%% Scene Alone (1s default)
            Screen('DrawTexture',dv.disp.ptr,dv.trial.sceneImageTexture); % Draw Texture (same texture variable for both trial types?)
            
            if dv.trial.ttime >= dv.pa.sceneTime + dv.trial.timeShowCueStart
                
                if strcmp(dv.trialType, 'study')
                    dv.trial.state = dv.states.SHOWPAIR;
                elseif strcmp(dv.trialType, 'test')
                    dv.trial.state = dv.states.DELAY;
                elseif dv.pa.strictDelay == 1 && strcmp(dv.trialType, 'test')
                    dv.trial.state = dv.states.FPONDELAY; % added strict delay
                end
            end
            
        end
        
    end

% SHOW PAIR
%---------------------------------------------------------------------%

    function dv = showPair(dv)
        if dv.trial.state == dv.states.SHOWPAIR 
            % Just for study trials
            
            dv.trial.fixFlagOn = 0;
            dv.trial.fixDotOn = 0;
            
            if dv.trial.showPairFlag
            dv.trial.ttime = GetSecs - dv.trial.trstart;
            dv.trial.timeShowPairStart = dv.trial.ttime;
            dv.trial.showPairFlag = 0;
            end
            
            %%%%%% Scene and Object (2s default)
            % Scene
            Screen('DrawTexture',dv.disp.ptr,dv.trial.sceneImageTexture);
            % Object
            Screen('DrawTexture',dv.disp.ptr,dv.trial.objectImageTexture,[],dv.trial.destRect); % Draw Texture
            
            if dv.trial.ttime >= dv.pa.showPairTime + dv.trial.timeShowPairStart
                dv.trial.ttime = GetSecs - dv.trial.trstart;
                dv.trial.timeComplete = dv.trial.ttime;
                dv.trial.state = dv.states.TRIALCOMPLETE;
            end
            
        end %end if trial state
    end % end func

 % DELAY 
 %---------------------------------------------------------------------%
    function dv = delay(dv)
        if dv.trial.state == dv.states.DELAY
            
            dv.trial.fixFlagOn = 1; 
            dv.trial.fixDotOn = 1;
            
            if dv.trial.delayFlag
                dv.trial.ttime = GetSecs - dv.trial.trstart;
                dv.trial.timeDelayStart = dv.trial.ttime;
                dv.trial.delayFlag = 0;
            end
            
            % no longer using box
%             if dv.trial.delayRectOn
%                 Screen('FillRect', dv.disp.ptr, dv.pa.delayBoxColor, dv.pa.centeredDelayRect);
%             end
            
            % need a grace time for them to look back at the center
            if dv.trial.ttime < dv.trial.graceTime + dv.trial.timeDelayStart
                %dv.trial.delayRectOn = 1;
                dv.trial.fixFlagOn = 1;
                dv.trial.fixDotOn = 1;
                
            elseif dv.trial.ttime < dv.trial.timeDelayStart + dv.trial.graceTime + dv.pa.delayTime && fixationHeld(dv)  % && dv.trial.ttime < dv.pa.delayTime + dv.trial.timeDelayStart + dv.trial.graceTime
                %%%% Delay (4.5s default)
                %dv.trial.delayRectOn = 1; % redundant (already set by this point)
                dv.trial.fixFlagOn = 1;
                dv.trial.fixDotOn = 1;
                
                % briefly represent fixation pt to cue onset of memory
                % probe (.5s)
            elseif dv.trial.ttime < dv.pa.delayTime + dv.trial.timeDelayStart + dv.trial.graceTime + dv.pa.probeCueTime && fixationHeld(dv)
                
                % Red dot off for visual cue
                dv.trial.fixDotOn = 0;
                dv.trial.fixFlagOn = 1;
                %dv.trial.delayRectOn = 0;
                
            elseif dv.trial.ttime >= dv.pa.delayTime + dv.trial.timeDelayStart + dv.pa.probeCueTime + dv.trial.graceTime && fixationHeld(dv)
                dv.trial.state = dv.states.SHOWPROBE;
                
            elseif ~fixationHeld(dv) && dv.trial.ttime > dv.pa.delayTime + dv.trial.timeDelayStart + dv.pa.probeCueTime + dv.trial.graceTime
                dv.trial.timeBreakFix = GetSecs - dv.trial.trstart; % always have to mark the time of break fixation
                dv.trial.state = dv.states.BREAKFIX; 
            end
        end
    end
    
% SHOW PROBE
%---------------------------------------------------------------------%
    function dv = showProbe(dv)
        if dv.trial.state == dv.states.SHOWPROBE
            % just for test trials
            
            dv.trial.fixFlagOn = 0;
            dv.trial.fixDotOn = 0;
            
            if dv.trial.showProbeFlag
                dv.trial.ttime = GetSecs - dv.trial.trstart;
                dv.trial.timeShowProbeStart = dv.trial.ttime;
                dv.trial.showProbeFlag = 0;
            end
            
            %%%%%% Scene and 2 Objects (2s)
            % ******** Scene ********
            Screen('DrawTexture',dv.disp.ptr,dv.trial.sceneImageTexture);
            
            % ******** Object 1 ********
            Screen('DrawTexture',dv.disp.ptr,dv.trial.object1ImageTexture,[],dv.trial.destRect1); % Draw Texture
            
            % ******** Object 2 ********
            Screen('DrawTexture',dv.disp.ptr,dv.trial.object2ImageTexture,[],dv.trial.destRect2); % Draw Texture
            
            if dv.trial.ttime >= dv.pa.probeTime + dv.trial.timeShowProbeStart
                dv.trial.ttime = GetSecs - dv.trial.trstart;
                dv.trial.timeComplete = dv.trial.ttime;
                dv.trial.state = dv.states.TRIALCOMPLETE;
            end
            
            % Quick and dirty eye measures
            dv.trial.eyePosi = dv.trial.eyePosi + 1;
            if dv.trial.eyeX >= dv.trial.destRect1(1) && dv.trial.eyeX <= dv.trial.destRect1(3) && dv.trial.eyeY >= dv.trial.destRect1(2) && dv.trial.eyeY <= dv.trial.destRect1(4) 

                % need to add in separate timer and easy first saccade
                % metric (don't have to though can calcluate everything)
%             dv.trial.timeRect1Start = GetSecs;
%             dv.trial.cumTimeRect1 = GetSecs - dv.trial.timeRect1Start;
%             
            dv.trial.eyeLoc(dv.trial.eyePosi) = 1;
            
            elseif dv.trial.eyeX >= dv.trial.destRect2(1) && dv.trial.eyeX <= dv.trial.destRect2(3) && dv.trial.eyeY >= dv.trial.destRect2(2) && dv.trial.eyeY <= dv.trial.destRect2(4)
            dv.trial.eyeLoc(dv.trial.eyePosi) = 2;
            
            else
            dv.trial.eyeLoc(dv.trial.eyePosi) = 0;    
            end
            
        end
    end





end % end master function