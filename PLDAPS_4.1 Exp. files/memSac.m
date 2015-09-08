function p = memSac(p,state)
% trial file for memSac (memory guided saccades) 
% relies on memSacSetup.m which sets up initial parameters

% Runs through default trial states (Drawing, frame flipping and updating etc.)
pldapsDefaultTrialFunction(p,state);

% *** TRIAL STATES *** Only called once at the beginning of each trial                                     
 switch state
    case p.trial.pldaps.trialStates.trialSetup
        % setup task
        %p.trial.stimulus.task = (rand() < p.trial.stimulus.proportionMemory) + 1; % Only testing memory guided right now
        if p.trial.stimulus.task == 1
            p.trial.stimulus.winScale = p.trial.stimulus.winScaleVisual;
        else
            p.trial.stimulus.winScale = p.trial.stimulus.winScaleMemory;
        end
        
        % ---Setup Targets--- %
        % Calculate Targets locations
        p.trial.stimulus.targ1XY = [p.trial.stimulus.targ1XYdeg(1)*p.trial.display.w2px(1) p.trial.stimulus.targ1XYdeg(2)*p.trial.display.w2px(2)]; 
        p.trial.stimulus.targ2XY = [p.trial.stimulus.targ2XYdeg(1)*p.trial.display.w2px(1) p.trial.stimulus.targ2XYdeg(2)*p.trial.display.w2px(2)];
        % Randomize which target is shown
        if rand() > .5     
        p.trial.stimulus.targLoc = p.trial.stimulus.targ1XY;
        else    
        p.trial.stimulus.targLoc = p.trial.stimulus.targ2XY;
        end
        % Calculate window of target
        p.trial.stimulus.targWin = p.trial.stimulus.fpWin + [1 1] * (sqrt(sum((p.trial.stimulus.targLoc - p.trial.stimulus.fixationXY).^2))*p.trial.stimulus.winScale);
        p.trial.stimulus.targRect = [-p.trial.stimulus.targWin(1) -p.trial.stimulus.targWin(2) p.trial.stimulus.targWin(1) p.trial.stimulus.targWin(2)] + [p.trial.display.ctr(1:2) + p.trial.stimulus.targLoc,p.trial.display.ctr(1:2) + p.trial.stimulus.targLoc];
        %--------------
        
        % Starting stimulus states for fixaiton (p.trial.stimulus.states
        % different then p.trial.pldaps.trialStates)
        p.trial.state = p.trial.stimulus.states.START; % NEED THIS!
        p.trial.stimulus.showFixationPoint = 1; % fix on (could just put this in setup)
        p.trial.stimulus.timeFpOn = NaN; % do we have to explicitly NaN these???? i don't think so
        p.trial.stimulus.timeTargetOn = NaN;
        %p.trial.stimulus.showTarg = 0; %reset targ, not needed now
        
        p.trial.stimulus.fpWin = p.trial.stimulus.fpWin*p.trial.display.ppd; % deg to pix
        
        
        % *** FRAME STATES *** % 
    case p.trial.pldaps.trialStates.framePrepareDrawing;
         
        p = checkFixation(p);
        
        p = checkTargetFixation(p);
        
        p = checkTrialState(p);
        
        % save eye position on every frame
        p.trial.stimulus.eyeXYs(1:2,p.trial.iFrame)= [p.trial.eyeX-p.trial.display.pWidth/2; p.trial.eyeY-p.trial.display.pHeight/2];
        
     case p.trial.pldaps.trialStates.frameDraw;

         % Drawing Fixation pt
         if p.trial.stimulus.showFixationPoint && p.trial.ttime > p.trial.stimulus.preTrial
             Screen('Drawdots',  p.trial.display.overlayptr,  p.trial.stimulus.fixationXY, ...
                 p.trial.stimulus.fixdotW , p.trial.display.clut.fixation, p.trial.display.ctr(1:2),1);
             if isnan(p.trial.stimulus.timeFpOn) % not sure if these time stamps are working
                 p.trial.stimulus.timeFpOn = p.trial.ttime;
                 p.trial.stimulus.frameFpOn = p.trial.iFrame;
             end
         end
         
         % Drawing Target window
         if p.trial.stimulus.showTargWin
             Screen('FrameRect', p.trial.display.overlayptr, p.trial.display.clut.window, p.trial.stimulus.targRect); % perhaps should be regular screen(not overlay) if actually showing it (which screen is which screen?)
         end
         
         % Drawing target 
         if ~isnan(p.trial.stimulus.timeFpEntered) %&& p.trial.stimulus.showTarg % flag not needed now
             if p.trial.ttime > (p.trial.stimulus.targOnset + p.trial.stimulus.timeFpEntered) && p.trial.ttime < (p.trial.stimulus.targOnset + p.trial.stimulus.timeFpEntered + p.trial.stimulus.targDuration(2)) % targ duration 2 is for memory guided, change this when vis-guided built in
                 Screen('Drawdots',  p.trial.display.overlayptr, p.trial.stimulus.targLoc, ...
                     p.trial.stimulus.fixdotW, p.trial.display.clut.red, p.trial.display.ctr(1:2),1);
                 if isnan(p.trial.stimulus.timeTargetOn) % not sure if these time stamps are working
                     p.trial.stimulus.timeTargetOn = p.trial.ttime;
                     p.trial.stimulus.frameTargetOn = p.trial.iFrame;
                 end
             end
         end
         
         %__________________________%
         
        % If you want the trial to have a maximum time based on frames
%     case p.trial.pldaps.trialStates.frameFlip;
%         % if the trial has exceeded maximum duration, end trial and move to
%         % the next
%         if p.trial.iFrame == p.trial.pldaps.maxFrames
%             p.trial.flagNextTrial=true;
%         end
        
         
 end

 
end % end memSac
 
%%helper functions
%-------------------------------------------------------------------------%
%%% INLINE FUNCTIONS
%-------------------------------------------------------------------------%
    function held = fixationHeld(p)
                held = squarewindow(p.trial.pldaps.pass,p.trial.display.ctr(1:2)+p.trial.stimulus.fixationXY-[p.trial.eyeX p.trial.eyeY],p.trial.stimulus.fpWin(1),p.trial.stimulus.fpWin(2));
    end
   
    function held = targetHeld(p)
                held = squarewindow(p.trial.pldaps.pass,p.trial.display.ctr(1:2)+p.trial.stimulus.targLoc-[p.trial.eyeX p.trial.eyeY],p.trial.stimulus.targWin(1),p.trial.stimulus.targWin(2));
    end
    
  
   
    
% CHECK FIXATION
%---------------------------------------------------------------------%
function p = checkFixation(p)
        % WAITING FOR SUBJECT FIXATION
        fixating=fixationHeld(p);
        if  p.trial.state == p.trial.stimulus.states.START
            
            if fixating && p.trial.ttime  < (p.trial.stimulus.preTrial+p.trial.stimulus.fixWait) && p.trial.ttime > p.trial.stimulus.preTrial % last part - to make sure the target is not accidently shown before the fp 
                
                p.trial.stimulus.timeFpEntered = p.trial.ttime;
                p.trial.stimulus.frameFpEntered = p.trial.iFrame;
                
                %p.trial.stimulus.showTarg = 1; % FLAG TARG, not needed anymore
                
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.FIXATION,p.trial.pldaps.iTrial)
                end
                p.trial.state = p.trial.stimulus.states.FPHOLD;
               
            elseif p.trial.ttime  > (p.trial.stimulus.preTrial+p.trial.stimulus.fixWait) 
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.BREAKFIX,p.trial.pldaps.iTrial)
                end
                p.trial.stimulus.timeBreakFix = p.trial.ttime;
                p.trial.state = p.trial.stimulus.states.BREAKFIX;
            end
        end
        
        % check if fixation is held
        
        if p.trial.state == p.trial.stimulus.states.FPHOLD
            
            if fixating && p.trial.ttime > (p.trial.stimulus.timeFpEntered + p.trial.stimulus.fpHoldTime) % || p.trial.iFrame==p.trial.pldaps.maxFrames) % still need max frames thing here?
                
                % Turn off fixation point (queue to make saccade)
                p.trial.stimulus.showFixationPoint = 0;
                
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.FIXATION,p.trial.pldaps.iTrial)
                end
%                 p.trial.ttime      = GetSecs - p.trial.trstart;
                p.trial.stimulus.timeFpOff  = p.trial.ttime;
                p.trial.stimulus.frameFpOff = p.trial.iFrame;
                
                p.trial.state = p.trial.stimulus.states.CHOOSETARG; % NEXT STATE: TARGET (need some grace time to saccade) 
                
            elseif ~fixating && p.trial.ttime < (p.trial.stimulus.timeFpEntered + p.trial.stimulus.fpHoldTime)                
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.BREAKFIX,p.trial.pldaps.iTrial)
                end
                p.trial.stimulus.timeBreakFix = GetSecs - p.trial.trstart;
                p.trial.state = p.trial.stimulus.states.BREAKFIX;
           
            end
        end
                
end    


% CHECK TARGET FIXATION
function p = checkTargetFixation(p)
    
    fixatingTarget = targetHeld(p);
    
    % CHOOSE TARGET
    if p.trial.state == p.trial.stimulus.states.CHOOSETARG
        
        if fixatingTarget % && p.trial.ttime % don't think I need anything else here
            p.trial.stimulus.timeTargEntered = p.trial.ttime;
            p.trial.state = p.trial.stimulus.states.HOLDTARG;
        elseif ~fixatingTarget && p.trial.ttime > (p.trial.stimulus.timeFpOff + p.trial.stimulus.targWait)
            p.trial.state = p.trial.stimulus.states.BREAKFIX;
        end
    end

    % HOLD TARGET
    if p.trial.state == p.trial.stimulus.states.HOLDTARG
        if fixatingTarget && p.trial.ttime > (p.trial.stimulus.timeTargEntered + p.trial.stimulus.targHold)
            % TRIALCOMPLETE
            p.trial.stimulus.timeComplete = p.trial.ttime; 
            p.trial.state = p.trial.stimulus.states.TRIALCOMPLETE;       
        elseif ~fixatingTarget && p.trial.ttime < (p.trial.stimulus.timeTargEntered + p.trial.stimulus.targHold)
            % BREAKFIX
            p.trial.state = p.trial.stimulus.states.BREAKFIX;
        end
    end
    
end

% TRIAL COMPLETE? -- GIVE REWARD IF GOOD
%---------------------------------------------------------------------%
function p = checkTrialState(p)
        if p.trial.state == p.trial.stimulus.states.TRIALCOMPLETE
            p.trial.pldaps.goodtrial = 1;
            
%             if p.trial.ttime > p.trial.stimulus.timeComplete + .2
                
              pds.behavior.reward.give(p);
              
%                 if p.trial.sound.use
%                     if p.trial.datapixx.use
%                         pds.datapixx.flipBit(p.trial.event.REWARD,p.trial.pldaps.iTrial);
%                     end
%                     PsychPortAudio('Start', p.trial.sound.reward, 1, [], [], GetSecs + .1);
%                     p.trial.ttime = GetSecs - p.trial.trstart;
%                     p.trial.stimulus.timeReward(p.trial.stimulus.iReward) = p.trial.ttime;
%                     p.trial.stimulus.iReward = p.trial.stimulus.iReward + 1;
%                     if p.trial.pldaps.pass == 0 && p.trial.datapixx.use
%                         pds.datapixx.analogOut(p.trial.stimulus.rewardTime)
%                     end
%                 end
%             end
            
%             if p.trial.ttime > p.trial.stimulus.timeComplete + 1
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.TRIALEND,p.trial.pldaps.iTrial);
                end
                p.trial.flagNextTrial = true;
%             end
            
        end
        
        if p.trial.state == p.trial.stimulus.states.BREAKFIX

            p.trial.pldaps.goodtrial = 0;

            if p.trial.sound.use && ~isnan(p.trial.stimulus.timeFpEntered) 
                PsychPortAudio('Start', p.trial.sound.breakfix, 1, [], [], GetSecs + .1);
            end
            
%             if p.trial.ttime > p.trial.stimulus.timeBreakFix + p.trial.stimulus.breakFixPenalty
                p.trial.flagNextTrial = true;
%             end
        end
        
        
end
    
    