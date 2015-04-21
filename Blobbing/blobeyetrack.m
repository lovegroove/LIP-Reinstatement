function dv = blobeyetrack(dv)

% **** uses normal runPLDAPS not cv version

% condition file
dv.trialFunction = 'runBlobEyeTrack';
dv.pass = 0;
dv.useMouse = 0; %also probably not needed
dv.showMouse = 0; % probably not even needed
dv.finish = 50; % # of trials
dv.trialTime = 10;      % trial duration in s
dv.randSig = 1;
dv.randStep = 1;
dv.choiceOn = 1; % flag for direction choice (did I implement this yet?)

dv = pdsDefaultTrialSetupEric(dv);

 %% Params
 
 % blob
 dv.thisSig = 8; % 8,10,13,16,19,22
 dv.sigRang = [8 10 13 16 19 22];
 dv.gabSize = 200;
 dv.cont = .5;
 dv.incr = 1;     % i^2 for dark, -i^2 for light
 dv.halfg = round(dv.gabSize/2);
 dv.stepRang = 300:100:800;
 dv.stepSize = 500;  % (1),(10) ave or max stepsize in pixels - determines target speed ( / dv.disp.ppd for deg/s)

 
 % noise
 dv.noysSize = 600; %300
 dv.walkScale = 200; % (100) effectively the half-width of the bounding box of the blob when doing velocity walks (trying 200)
 dv.noysCont = 1;  % noys contrast
 
 % states
 dv.states.BLOBWALK = 1;
 dv.states.CHOICE = 2;
 













end
