function dv = blobeyetrack(dv)

% **** uses normal runPLDAPS not cv version

% condition file
dv.trialFunction = 'runBlobEyeTrack';
dv.pass = 0;
dv.useMouse = 0;
dv.showMouse = 0; % probably not even needed
dv.finish = 2; %1e3; % # of trials
dv.trialTime = 8;      % trial duration in s
dv.choiceOn = 1; % flag for direction choice (did I implement this yet?)

dv = pdsDefaultTrialSetupEric(dv);

 %% Params
 
 % blob
 dv.thisSig = 10; % 8,10,13,16,19,22
 dv.gabSize = 200;
 dv.cont = .5;
 dv.incr = 1;     % i^2 for dark, -i^2 for light
 dv.halfg = round(dv.gabSize/2);
 dv.stepSize = 10;      % (1) ave or max stepsize in pixels - determines target speed

% noise
dv.noysSize = 500; %300
dv.walkScale = 100; % effectively the half-width of the bounding box of the blob when doing velocity walks
dv.noysCont = 1;  % noys contrast

% states
dv.states.BLOBWALK = 1; 
dv.states.CHOICE = 2;














end
