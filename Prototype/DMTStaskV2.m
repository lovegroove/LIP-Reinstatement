% Memory Task in the PPC - can set the trial type to "study" or "test"
% VERSION 2
% Relies upon functions: loadImagePairsFunc and calcVisAng
clear all;
close all;
sca;

% load images function
[pairOrder, fileInfo] = loadImagePairsFunc;

%% Psychtoolbox startup stuff

% Call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2); %Only use if you are using the PsychImaging function
Priority(1);

% Unify Keys
KbName('UnifyKeyNames');

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen
screenNumber = max(screens);

% Remove startup stuff
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
inc = white - grey;

% Open a screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey); %Pscyh Imaging func or just Screen

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Get the center coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Fixation Cross and fixation rectangle and delay box

% Setup the text type for the window
Screen('TextFont', window, 'Ariel');
Screen('TextSize', window, 36);

% Set the size of the arms of our fixation cross
fixCrossDimPix = 40;

% Set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];

% Set the line width for our fixation cross
lineWidthPix = 4;

% Make a base Rect of 200 by 200 pixels
baseRect = [0 0 200 200];

% Make Delay box
delayRect = [0 0 50 50];
centeredDelayRect = CenterRectOnPointd(delayRect, xCenter, yCenter);

% Fixation Rectangle
fixationRect = [0 0 100 100];
centeredRect = CenterRectOnPointd(fixationRect, xCenter, yCenter);

%%  Mouse cursor & Trial while loop Parameters
% Initialize stuff
stopkey = KbName('ESCAPE');
pausekey = KbName('p');
trialIndex = 1;
blockIndex = 1;
trialCumulativeTimer = 0;
i = 1; % Just for incremental testing
fixFlag = 0;
virgin = 1;
idleTime = 0;

%%%%% Trial Parameters and Timing *********** THINGS YOU MIGHT WANT TO CHANGE ***********
trialType = 'test'; % study or test
numTrials = 42;  %numTrials = length(objectImages_mask); % normally this?
correctObjectLocs = cell(1,numTrials); % For saving correct objects and their placement
showMouse = 0;
showFixWindow = 0;
useAudio = 0;
showLagReport = 0;
cursorColor = black;
fixWindowColor = white;
fixCrossColor = white;
objectSize = .25; % In percent of base rectangle

% Sub-trial Event Timing (in ms) (if random timing prefered, check virgin if statement within loop)
fixPtTime = 2000; % currently randi(5) in loop (temporarily testing no variability)
sceneTime = 1000;
objectSceneTime = 2000;
probeCueTime = 500;
delayTime = 5000;
graceTime = 500;

%%%%% Test Trial Stimulus Location and Geometry Parameters
alpha = 6; % degrees of visual angle
Z = 57; % in cm
scrX = 16; % in inches
scrY = 9; % in inches
[Dx, Dy] = calcVisAng(alpha, Z, scrX, scrY, screenXpixels, screenYpixels); % Func to calculate visual angle

% Stimulus (object) angle and distance from center (fixation pt)
% calculations (can vary this randomly or systemactically if desired (put in and opposite RF), *** must put in the loop in the random case***
theta = 135 * (pi/180); % in degrees, convert to radians
r = sqrt((Dx)^2 + (Dy)^2); % degrees of visual angle (which r will be the x and y dimensions)
object1loc = [xCenter + r * cos(theta) yCenter + r * sin(theta)];
object2loc = [xCenter + r * cos(theta+pi) yCenter + r * sin(theta+pi)];

%%%%% End Test Trial Params

%% Audio Setup
beepBreakFix = 'C:\Users\MrLovegroove\Documents\MATLAB\beepsounds\breakfix';
beepTrialStart = 'C:\Users\MrLovegroove\Documents\MATLAB\beepsounds\cue';

InitializePsychSound;

%% Main Loop

% Deal with cursor & snyc flip
%SetMouse(xCenter, yCenter, window);
HideCursor;
vbl = Screen('Flip', window); %Initially synchronize with retrace, take base time in vbl

%%%%%%%%% Trial Loop%%%%%%%%%%%%%%%%
while trialIndex <= numTrials % make sure you remember this might depend on trial type*******
    
    startTime = GetSecs;
    
    % Shuffle order every 40 trials (each time they see all of the pairings
    % once) and add more
    if blockIndex == 40
        pairOrder = [pairOrder; datasample(pairOrder,40,1,'replace',false)]; % Does this work ok????
        blockIndex = 0;
    end
    
    % Menu Options
    [keyIsDown,secs,keyCode] = KbCheck;
    % ESC was pressed stop display
    if keyCode(stopkey)
        disp('ESC pressed, exiting trial.');
        break;
    elseif keyCode(pausekey)
        keyboard;
        disp('Trial paused, Type RETURN to continue execution.');
    end
    
    % Get mouse position
    [x,y,buttons] = GetMouse(window);
    
    % See if the mouse cursor is inside the fixation window
    insideFix = IsInRect(x, y, centeredRect);
    
    switch trialType     % number of trials might depend on type as well*****
                %%%%%%%%%%%%%
                %Study Block%**********************************************
                %%%%%%%%%%%%%
        case 'study'
            % Trial execution block (1 loop =~16.6 ms (frame rate 60hz?))
            if virgin %Front load computationally heavy shit, so minor lag is only between trials
                virgin = 0;
                startStudyTrial = GetSecs;
                % Set random durations for events if perferred 
                %fixPtTime = 1000 * randi(5); %also acts as ITI *********** No random interval now for testing
                
                % Load in Image files
                objectImageFile = [fileInfo.objectPath pairOrder{trialIndex,1}];
                sceneImageFile = [fileInfo.scenePath pairOrder{trialIndex,2}];  
                
                % Background
                sceneImage = imread(sceneImageFile, fileInfo.sceneFormat); % Read in image
                sceneImage = imresize(sceneImage, [screenYpixels screenXpixels]); % Resize Image to screen
                sceneImageTexture = Screen('MakeTexture',window,sceneImage); % Create Texture
                
                % Object
                objectImage = imread(objectImageFile, fileInfo.objectFormat); % Read in image
                objectImageTexture = Screen('MakeTexture',window,objectImage); % Create Texture
                [s1, s2, s3] = size(objectImage);
                baseRect = [0 0 s1 s2];
                destRect = CenterRectOnPointd(baseRect .* objectSize, xCenter, yCenter); %Set size and location

                % Beep Audio - Break Fixation
                if useAudio
                    [wavdata, freq] = wavread(beepTrialStart);
                    wavdata = wavdata';
                    nrchannels = size(wavdata,1);
                    pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
                    PsychPortAudio('FillBuffer', pahandle, wavdata);
                    PsychPortAudio('Start', pahandle, [], 0, 1);
                    %PsychPortAudio('Stop', pahandle);
                end
                
                % Wait buffer between trials
                WaitSecs(.5);
                
            elseif trialCumulativeTimer <= fixPtTime  % need to implement while graceTime < fixationTime
                %%%%%% Fixation cross (1-5s) testing 2s
                Screen('DrawLines', window, allCoords, lineWidthPix, fixCrossColor, [xCenter yCenter], 2); % Draw the fixation Cross
                if insideFix && ~fixFlag % First time in FW
                    %beep;
                    fixFlag = 1;
                elseif ~insideFix && fixFlag && trialCumulativeTimer < graceTime
                    % Continue (As you were soldier!)
                elseif ~insideFix && fixFlag && trialCumulativeTimer > graceTime % Left FW and over grace time, how to use GetSecs here?
                    % Beep Audio - Break Fixation
                    if useAudio
                        [wavdata, freq] = wavread(beepBreakFix);
                        wavdata = wavdata';
                        nrchannels = size(wavdata,1);
                        pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
                        PsychPortAudio('FillBuffer', pahandle, wavdata);
                        PsychPortAudio('Start', pahandle, [], 0, 1);
                        %PsychPortAudio('Stop', pahandle);
                    end
                    %SetMouse(xCenter, yCenter, window);
                    %WaitSecs(.5);
                    disp('Broke Fixation');
                    fixFlag = 0;
                    trialCumulativeTimer = 0;
                elseif ~insideFix && ~fixFlag % if he/she is never looking at it
                    trialCumulativeTimer = 0;
                    fixFlag = 0;
                    idleTime = idleTime + loopTime * 1000;
                end
            elseif trialCumulativeTimer > fixPtTime && trialCumulativeTimer <= sceneTime + fixPtTime && fixFlag % and fixFlag here
                %%%%%% Scene alone (1s)
                Screen('DrawTexture',window,sceneImageTexture); % Draw Texture
                
            elseif trialCumulativeTimer > sceneTime + fixPtTime && trialCumulativeTimer <= fixPtTime + sceneTime + objectSceneTime
                %%%%%% Scene and Object (2s)
                % Scene
                Screen('DrawTexture',window,sceneImageTexture);
                % Object
                Screen('DrawTexture',window,objectImageTexture,[],destRect); % Draw Texture
                
            elseif trialCumulativeTimer > fixPtTime + sceneTime + objectSceneTime
                % Increment trial and reset timer
                trialIndex = trialIndex + 1;
                blockIndex = blockIndex + 1;
                trialCumulativeTimer = 0;
                virgin = 1;
                endStudyTrial = GetSecs - startStudyTrial;
                %fixPtTime = randi(5) * 1000 / loopTime; %why doesn't this work
                Screen('Close');
                %WaitSecs(.5);
            end %end study trial block
            
        case 'test' 
                %%%%%%%%%%%%
                %Test Block%***********************************************
                %%%%%%%%%%%%
            if virgin 
                virgin = 0;
                startTestTrial = GetSecs;
                % Set random durations for events if perferred 
                %fixPtTime = 1000 * randi(5);
                
                % Load in Image files - Synced up with Test Trial***** change objectImageFile to correct...
                objectImageFile = [fileInfo.objectPath pairOrder{trialIndex,1}];
                sceneImageFile = [fileInfo.scenePath pairOrder{trialIndex,2}];
                
                % Randomly select foil and place object images randomly - can choose number of objects here
                foilObject = pairOrder{randsample(length(pairOrder),1) ,1};
                
                % Also reroll here just in case you happen to pick the same object
                while strcmp(foilObject, pairOrder{trialIndex,1})
                    foilObject = pairOrder{randsample(length(pairOrder),1) ,1};
                end
                % Place objects randomly
                foilObjectImageFile = [fileInfo.objectPath foilObject];
                placedObjects = cell(2,1);
                placedObjects{1} = objectImageFile;
                placedObjects{2} = foilObjectImageFile;
                [object1 object2] = placedObjects{randsample(2,2)};
                
                % Background - make texture
                sceneImage = imread(sceneImageFile, fileInfo.sceneFormat); % Read in image
                sceneImage = imresize(sceneImage, [screenYpixels screenXpixels]); % Resize Image to screen
                sceneImageTexture = Screen('MakeTexture',window,sceneImage); % Create Texture
                
                % Object 1 - make texture and assign location
                object1Image = imread(object1, fileInfo.objectFormat); % Read in image
                object1ImageTexture = Screen('MakeTexture',window,object1Image); % Create Texture
                [s1, s2, s3] = size(object1Image);
                baseRect = [0 0 s1 s2];
                destRect1 = CenterRectOnPointd(baseRect .* objectSize, object1loc(1), object1loc(2)); %Set size and location
                
                % Object 2 - make texture and assign location
                object2Image = imread(object2, fileInfo.objectFormat); % Read in image
                object2ImageTexture = Screen('MakeTexture',window,object2Image); % Create Texture
                [s1, s2, s3] = size(object2Image);
                baseRect = [0 0 s1 s2];
                destRect2 = CenterRectOnPointd(baseRect .* objectSize, object2loc(1), object2loc(2)); %Set size and location
                
                % Save correct object and location 
                correctObjectLocs{trialIndex} = {placedObjects{1}, placedObjects{2}; destRect1, destRect2};
                
                % Beep Audio - Break Fixation
                if useAudio
                    [wavdata, freq] = wavread(beepTrialStart);
                    wavdata = wavdata';
                    nrchannels = size(wavdata,1);
                    pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
                    PsychPortAudio('FillBuffer', pahandle, wavdata);
                    PsychPortAudio('Start', pahandle, [], 0, 1);
                    %PsychPortAudio('Stop', pahandle);
                end
                
                % Wait buffer between trials
                WaitSecs(.5);
                
            elseif trialCumulativeTimer <= fixPtTime % Fixation Point
                %%%%%% Fixation cross (1-5s) testing 2s
                Screen('DrawLines', window, allCoords, lineWidthPix, fixCrossColor, [xCenter yCenter], 2); % Draw the fixation Cross
                 if insideFix && ~fixFlag % First time in FW
                    %beep;
                    fixFlag = 1;
                elseif ~insideFix && fixFlag && trialCumulativeTimer < graceTime
                    % Continue (As you were soldier!)
                elseif ~insideFix && fixFlag && trialCumulativeTimer > graceTime % Left FW and over grace time, how to use GetSecs here?
                    % Beep Audio - Break Fixation
                    if useAudio
                        [wavdata, freq] = wavread(beepBreakFix);
                        wavdata = wavdata';
                        nrchannels = size(wavdata,1);
                        pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
                        PsychPortAudio('FillBuffer', pahandle, wavdata);
                        PsychPortAudio('Start', pahandle, [], 0, 1);
                        %PsychPortAudio('Stop', pahandle);
                    end
                    %SetMouse(xCenter, yCenter, window);
                    %WaitSecs(.5);
                    disp('Broke Fixation');
                    fixFlag = 0;
                    trialCumulativeTimer = 0;
                elseif ~insideFix && ~fixFlag % if he/she is never looking at it
                    trialCumulativeTimer = 0;
                    fixFlag = 0;
                    idleTime = idleTime + loopTime * 1000;
                end
            elseif trialCumulativeTimer > fixPtTime && trialCumulativeTimer <= sceneTime + fixPtTime && fixFlag
                %%%%%% Scene alone (1s)
                Screen('DrawTexture',window,sceneImageTexture); % Draw Texture
                
            elseif trialCumulativeTimer > sceneTime + fixPtTime && trialCumulativeTimer <= delayTime + fixPtTime + sceneTime
                %%%% Delay (6s)
                Screen('FillRect', window, black, centeredDelayRect);
                
            elseif trialCumulativeTimer > delayTime + fixPtTime + sceneTime && trialCumulativeTimer <= delayTime + fixPtTime + sceneTime + probeCueTime
                % Switch back to fixation cross for visual cue of upcoming memory probe (.5s)
                Screen('DrawLines', window, allCoords, lineWidthPix, fixCrossColor, [xCenter yCenter], 2);
                
            elseif trialCumulativeTimer > delayTime + fixPtTime + sceneTime + probeCueTime && trialCumulativeTimer <= delayTime + fixPtTime + sceneTime + probeCueTime + objectSceneTime
                %%%%%% Scene and 2 Objects (2s)
                % Scene
                Screen('DrawTexture',window,sceneImageTexture);
                
                % ******** Object 1 ********
                Screen('DrawTexture',window,object1ImageTexture,[],destRect1); % Draw Texture
                
                % ******** Object 2 ********
                Screen('DrawTexture',window,object2ImageTexture,[],destRect2); % Draw Texture
                
            elseif trialCumulativeTimer > delayTime + fixPtTime + sceneTime + probeCueTime + objectSceneTime
                % Increment trial and reset timer
                trialIndex = trialIndex + 1;
                blockIndex = blockIndex + 1;
                trialCumulativeTimer = 0;
                virgin = 1;
                endTestTrial = GetSecs - startTestTrial;
                Screen('Close');
                %WaitSecs(.5);
            end % end ifelse test block
    end % end case - trial block
    
    
    %*********** Show Options ************
    % Draw a white dot where the mouse is
    if showMouse
        Screen('DrawDots', window, [x y], 10, cursorColor, [], 2);
    end
    % Draw the rect to the screen
    if showFixWindow
        Screen('FrameRect', window, fixWindowColor, centeredRect);
    end
    
    % FLIP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    vbl = Screen('Flip', window, vbl + 0.5*ifi);
    
    % Increment timer
    %WaitSecs(.001);
    loopTime = GetSecs - startTime;
    T(i) = loopTime;
    trialCumulativeTimer = trialCumulativeTimer + loopTime * 1000;
    i = i + 1;
end % end mouse loop

if showLagReport
% Checking latency issues
plot(T)
title('Latency (.5s of spikes is deliberate buffer ITI time)')
xlabel('Loop Iterations')
ylabel('Time (s)')
end

%FlushEvents;
ShowCursor;

if useAudio
% Stop and close the audio device
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
end

% Return PTB settings and clean up
Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
Screen('CloseAll')
sca;
Priority(0);

