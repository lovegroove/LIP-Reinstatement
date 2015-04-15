function [PDS,dv] = runBlobEyeTrack(PDS,dv)
% trial function for blob tracking with eyes

dv = defaultTrialVariables(dv); % setup default trial struct

%-------------------------------------------------------------------------%

% trial variables - MUST UPDATE ALL VARS
 dv.frameRate = round(1./Screen('GetFlipInterval', dv.disp.ptr));
    dv.totFrames = dv.trialTime*dv.frameRate;
    
    % make a Gaussian matlab matrix
    dv.alphGab = ggaus(dv.gabSize, dv.thisSig);   
    dv.targGab = dv.alphGab./(2.*pi.*dv.thisSig.^2);  % keep total energy constant
    dv.targGab = dv.cont.*2.*pi.*64*dv.targGab;               % scale a gaus with a sig of 8 has height 1
    dv.gabTex = Screen('MakeTexture', dv.disp.ptr, dv.targGab, [], [], 2);  % store as 32bit texture

% making the noys
dv.noysNormFac = dv.noysCont*(1/12);
    dv.theNoys= randn(dv.noysSize);
    dv.theNoys(abs(dv.theNoys)>3)= 3;
    dv.theNoys = dv.theNoys*dv.noysNormFac +0.5; % should put us in 0.25 to 0.75 range
    dv.noysTex = Screen('MakeTexture', dv.disp.ptr, dv.theNoys, [], [], 2);  % store as 32bit texture

  
    % random walk in velocity
    dv.walkCounter = 0; %initialize
    dv.walkState = rng;            
    temp = randn(2,dv.totFrames);
    dv.targCoords = dv.stepSize*cumsum(temp, 2);
    dv.targCoords = detrend(dv.targCoords');     
    dv.targCoords = cumsum(dv.targCoords', 2);  % positions
    dv.targCoords = dv.walkScale.*dv.targCoords./mmax(dv.targCoords);
    
    % WHAT ABOUT VARS FROM HERE DOWN
    
    % Move the cursor to the center of the screen
%     centerX = round(theRect(RectRight)/2);
%     centerY = round(theRect(RectBottom)/2);
%     dv.targCoords(1,:) = dv.targCoords(1,:) + centerX;
%     dv.targCoords(2,:) = dv.targCoords(2,:) + centerY;

    % will this not be the center of the rect no and be the center of the
    % screen? yes, does it matter?
    dv.targCoords(1,:) = dv.targCoords(1,:) + dv.pa.xCenter;
    dv.targCoords(2,:) = dv.targCoords(2,:) + dv.pa.yCenter;
    SetMouse(dv.pa.xCenter,dv.pa.yCenter);
    HideCursor();
    
    % Use key bindings from PLDAPS
%     % Keys
%     stopKey=KbName('esc');
%     leftKey=KbName('left');   
%     rightKey=KbName('right'); 

%-------------------------------------------------------------------------%

% % preallocate data aquisition variables
%flipTimes       = zeros(1e4,2);
loopTimes       = zeros(1e4,2);
%photodiodeTimes = zeros(1e4,2);
eyepos          = zeros(1e4,4);     % preallocate eye position data

% if dv.pass == 0
% dv.disp.inputSamplingRate = 240;
% eyeTimeStep     = 1/dv.disp.inputSamplingRate;
% end

%-------------------------------------------------------------------------%
%%% Degrees to pixels %%%
% We input positions on the screen in degrees of visual angle, but
% PsychToolbox takes calls to pixel locations. dv.disp.ppd is the pixels
% per degree of this monitor, given the viewing distance. It is important
% to set this up correctly in the setupPLADPS script and it requires
% Init_StereoDispPI to be implemented properly.

%dv.trial.fpWin = dv.pa.fpWin*dv.disp.ppd;  % commenting out for right now,
%do I need fpWin?

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

% starting state****
dv.trial.state = dv.states.BLOBWALK;

%% Main While Loop for the Trial******************************************
%---------------------------------------------------------------------%
while ~dv.trial.flagNextTrial && dv.quit == 0
    % trial time
    dv.trial.ttime  = GetSecs - dv.trial.trstart;
    
    dv.walkCounter = dv.walkCounter + 1;
    
    %%% get inputs %%%
    %---------------------------------------------------------------------%
    % get mouse/eyetracker data
    
    [dv.trial.cursorX,dv.trial.cursorY,dv.trial.isMouseButtonDown] = GetMouse;
    
    dv = pdsGetEyePosition(dv); % Eye position or Mouse Position if UseMouse is flagged (within func)
    %dv = pdsGetEyePosition(dv, updateQueue); % Need this update queue?????????????????
    
    %%% TRIAL STATES %%% dv.trial.state
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    % trial start states????
    
    
    % Blob tracking
    dv = blobWalk(dv);
    
    % Direction choice
    dv = dirChoice(dv);
    
    % trial complete state ????
    
    
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
    end
    
    
    %% DRAWING
    %---------------------------------------------------------------------%
    % Show Mouse pointer if desired
    if dv.showMouse
        Screen('DrawDots', dv.disp.ptr, [dv.trial.cursorX,dv.trial.cursorY], 10, [1,1,1], [], 2);
    end
    
    % FLIP
    vbl = Screen('Flip', dv.disp.ptr, vbl + 0.5*ifi);
    
    % needed here or after loop
    Screen('Close', dv.noysTex);  % issue dping this in loop*******????
    
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


 %-------------------------------------------------------------------------%
% PDS SAVING STUFF
if dv.pass == 0
    %Eyelink
    PDS.data.eyeposLoop{dv.j} = eyepos;
    PDS.data.eyelinkSampleBuffer{dv.j}  =  dv.el.sampleBuffer(:,~isnan(dv.el.sampleBuffer(1,:)));
    PDS.data.eyelinkEventBuffer{dv.j}   =  dv.el.eventBuffer(:,~isnan(dv.el.eventBuffer(1,:)));
    PDS.data.eyelinkQueueStruct = eye;
    
    PDS.eyepos{dv.j}         = eyepos(1:dv.trial.iSample+1,:);  % remove extra    
end

PDS.data.dirChoice{dv.j} = dv.dirChoice;



 %-------------------------------------------------------------------------%
 %%%%%%% INLINE FUNCTIONS: STATES
 %-------------------------------------------------------------------------%

    function dv = blobWalk(dv)
        if dv.trial.state == dv.states.BLOBWALK
            
            % it's blobbering time! - MUST UPDATE VARS!!!!
             % Draw NOYS!
            dv.theNoys= randn(dv.noysSize);
            dv.theNoys(abs(dv.theNoys)>3)= 3;
            dv.theNoys = dv.theNoys*dv.noysNormFac +0.5; % should put us in 0.25 to 0.75 range
            dv.noysTex = Screen('MakeTexture', dv.disp.ptr, dv.theNoys, [], [], 2);  % store as 32bit texture
            Screen('BlendFunction', dv.disp.ptr, GL_ONE, GL_ZERO);
            Screen('DrawTexture', dv.disp.ptr, dv.noysTex, [],[],[],[],0.5);
            
            
            % Draw Blob!
            Screen('BlendFunction', dv.disp.ptr, GL_SRC_ALPHA, GL_ONE);
            Screen('DrawTexture', dv.disp.ptr, dv.gabTex, [],...
                [dv.targCoords(1,dv.walkCounter), dv.targCoords(2,dv.walkCounter),...
                dv.targCoords(1,dv.walkCounter)+dv.gabSize, dv.targCoords(2,dv.walkCounter)+dv.gabSize]-dv.halfg, ...
                [],[],0.5);
            
            
            % end of trial switch to choice state
            if dv.walkCounter == dv.totFrames; % does this work?
            dv.trial.state = dv.states.CHOICE;
            end
        end
    end

    function dv = dirChoice(dv)
        if dv.trial.state == dv.states.CHOICE
            dv.trial.ttime = GetSecs - dv.trial.trstart; % probably need to record some times
            Screen(dv.disp.ptr,'DrawText','Is the blob in the right or left half of the box?  You have 3 seconds to respond.',50,50,255); % need good coords*****
                %[x,y,buttons] = GetMouse(whichScreen);
                if dv.trial.firstPressQ(dv.kb.Rarrow)
                    dv.dirChoice = 1; % these will have to be incremented when there are a block of trials
                    dv.trial.flagNextTrial = true;
                elseif dv.trial.firstPressQ(dv.kb.Larrow) %or instead dv.trial.firstPressQ(dv.kb.pKey) - whatever the name of the key is****************
                    dv.dirChoice = 2;
                    dv.trial.flagNextTrial = true;
%                 else
%                     direcKey = 0;
                end    
        end
    end


















end %end master function