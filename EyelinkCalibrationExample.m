% Short MATLAB example program that uses the Eyelink and Psychophysics
% Toolboxes to do calibration

clear all;
commandwindow;

if 1 Screen('Preference', 'SkipSyncTests', 1); end

fprintf('EyelinkToolbox Example\n\n\t');

% STEP 1
% Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.

if EyelinkInit()~= 1; %
    return;
end;
% STEP 2
% Open a graphics window on the main screen
% using the PsychToolbox's Screen function.

screenNumber=max(Screen('Screens'));
[window, wRect]=Screen('OpenWindow', screenNumber, 0,[],32,2);
Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% STEP 3
% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
el=EyelinkInitDefaults(window);

% make sure that we get gaze data from the Eyelink
Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');

[v vs]=Eyelink('GetTrackerVersion');

fprintf('Running experiment on a ''%s'' tracker.\n', vs );

% STEP 4
% Calibrate the eye tracker
EyelinkDoTrackerSetup(el);

% do a final check of calibration using driftcorrection
%EyelinkDoDriftCorrection(el); %For some reason this line prevents graceful exit from the program

