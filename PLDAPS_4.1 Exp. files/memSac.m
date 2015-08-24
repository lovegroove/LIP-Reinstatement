function p = memSac(p,state)
% trial file for memSac (memory guided saccades) 
% relies on memSacSetup.m which sets up initial parameters

% Runs through default trial states (Drawing, frame flipping and updating etc.)
pldapsDefaultTrialFunction(p,state);

% *** TRIAL STATES ***                                     
 switch state
    case p.trial.pldaps.trialStates.trialSetup
        % setup task
        %p.trial.stimulus.task = (rand() < p.trial.stimulus.proportionMemory) + 1; % Only testing memory guided right now
        if p.trial.stimulus.task == 1
            p.trial.stimulus.winScale = p.trial.stimulus.winScaleVisual;
        else
            p.trial.stimulus.winScale = p.trial.stimulus.winScaleMemory;
        end
        
        % setup Jonas' fancy framedrop plotting
        p.trial.sf = opticflow.screenFigure(p, [-p.trial.display.wWidth/2+10 -p.trial.display.wHeight/2+25], [20 5], p.trial.display.overlayptr, 12,[0 3000], [0 3*p.trial.display.ifi],p.trial.display.clut.hBlack*[1 1 1]);
        
        % ---Setup Targets---
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
        p.trial.stimulus.showFixationPoint = 1; % Or just do this in setup file...?
        
        
        % *** FRAME STATES *** 
    case p.trial.pldaps.trialStates.framePrepareDrawing;
        
        % MAKE SURE THIS IS WORKING, (remember next trial will be flagged
        % when max frames are reached) DOES ORDER OF THESE FUNCTIONS
        % MATTER??
        
        p = checkTrialState(p);
         
        p = checkFixation(p);
        
        p = checkTargetFixation(p);
        
        
     case p.trial.pldaps.trialStates.frameDraw;
         
         % Frame Drop plotting
         if p.trial.pldaps.draw.frameDropPlot.use
             %%draw stuff to the experimenter screen
             %we have to to it here, as, if we don't want it to draw over the dots drawn later
             dataInds=max(round(p.trial.iFrame-3/p.trial.display.ifi),1):p.trial.iFrame-1;
             time=dataInds*1000*p.trial.display.ifi;
             time=time(2:end);
             p.trial.sf.xlims=[max([min(time) 0]) max([time 3000])];
             opticflow.screenPlot(p.trial.sf, p.trial.sf.xlims, [p.trial.display.ifi p.trial.display.ifi], p.trial.display.clut.hBlack, '-');
             opticflow.screenPlot(p.trial.sf, time, diff(p.trial.timing.flipTimes(1,dataInds)), p.trial.display.clut.hRed, '--');
         end
         
%          %-------------------- if only switching colors-------------------%
%          % Drawing only (every frame, just switch colors)        
%          Screen('Drawdots',  p.trial.display.overlayptr,  p.trial.stimulus.fixationXY, ...
%              p.trial.stimulus.fixdotW , p.trial.stimulus.colorFixDot, p.trial.display.ctr(1:2),1);
%          
%          Screen('Drawdots',  p.trial.display.overlayptr, p.trial.stimulus.targ1XY, ...
%              p.trial.stimulus.fixdotW, p.trial.stimulus.colorTarget1Dot, p.trial.display.ctr(1:2),1);
%          
%          Screen('Drawdots',  p.trial.display.overlayptr, p.trial.stimulus.targ2XY, ...
%              p.trial.stimulus.fixdotW, p.trial.stimulus.colorTarget2Dot, p.trial.display.ctr(1:2),1);
         
         %__________________________%
         % Drawing Fixation pt
         if p.trial.stimulus.showFixationPoint % && p.trial.ttime > p.trial.stimulus.preTrial   
             Screen('Drawdots',  p.trial.display.overlayptr,  p.trial.stimulus.fixationXY, ...
                 p.trial.stimulus.fixdotW , p.trial.display.clut.fixation, p.trial.display.ctr(1:2),1);
         end
         
         % Drawing Target window
         if p.trial.stimulus.showTargWin
             Screen('FrameRect', p.trial.display.overlayptr, p.trial.display.clut.window, p.trial.stimulus.targRect); % perhaps should be regular screen(not overlay) if actually showing it
         end
         
         % Drawing target - Showing (flashing) the target does not have to
         % depend on target functions because subject is still fixating
         % when they are on
         if p.trial.stimulus.showTargets && p.trial.ttime > (p.trial.stimulus.targOnset + p.trial.stimulus.timeFpEntered) && p.trial.ttime < (p.trial.stimulus.targOnset + p.trial.stimulus.timeFpEntered + p.trial.stimulus.targDuration(2)) % targ duration 2 is for memory guided
             Screen('Drawdots',  p.trial.display.overlayptr, p.trial.stimulus.targLoc, ...
                 p.trial.stimulus.fixdotW, p.trial.display.clut.red, p.trial.display.ctr(1:2),1);    
         end
         
         % When to draw logic
%          if p.trial.ttime > (p.trial.stimulus.targOnset + p.trial.stimulus.preTrial) && ~p.trial.stimulus.showTargets
%              p.trial.stimulus.showTargets = 1;
%          end
%          
%          if p.trial.ttime > (p.trial.stimulus.targDuration(2) + p.trial.stimulus.targOnset + p.trial.stimulus.preTrial) && p.trial.stimulus.showTargets
%              p.trial.stimulus.showTargets = 0;
%          end
%          
%          if p.trial.ttime > (p.trial.stimulus.targWait + p.trial.stimulus.targDuration(2) + p.trial.stimulus.targOnset + p.trial.stimulus.preTrial)
%              p.trial.stimulus.showFixationPoint = 0;
%          end
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
    
    % Make sure this works right...
    function held = targetHeld(p)
                held = squarewindow(p.trial.pldaps.pass,p.trial.display.ctr(1:2)+p.trial.stimulus.targWin-[p.trial.eyeX p.trial.eyeY],p.trial.stimulus.targWin(1),p.trial.stimulus.targWin(2));
    end
    
  
   
    
% CHECK FIXATION
%---------------------------------------------------------------------%
function p = checkFixation(p)
        % WAITING FOR SUBJECT FIXATION
        fixating=fixationHeld(p);
        if  p.trial.state == p.trial.stimulus.states.START
            %p.trial.stimulus.fpWin = p.trial.fpWin*p.trial.display.ppd/2;
            if fixating && p.trial.ttime  < (p.trial.stimulus.preTrial+p.trial.stimulus.fixWait)
                
                p.trial.stimulus.timeFpEntered = p.trial.ttime;%GetSecs - p.trial.trstart;
                p.trial.stimulus.frameFpEntered = p.trial.iFrame;
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.FIXATION,p.trial.pldaps.iTrial)
                end
                p.trial.state = p.trial.stimulus.states.FPHOLD;
            elseif p.trial.ttime  > (p.trial.stimulus.preTrial+p.trial.stimulus.fixWait) 
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.BREAKFIX,p.trial.pldaps.iTrial)
                end
                p.trial.stimulus.timeBreakFix = p.trial.ttime;%GetSecs - p.trial.trstart;
                p.trial.state = p.trial.stimulus.states.BREAKFIX;
            end
        end
        
        % check if fixation is held
        
        if p.trial.state == p.trial.stimulus.states.FPHOLD
            %p.trial.stimulus.fpWin = p.trial.fpWin*p.trial.display.ppd;
            if fixating && p.trial.ttime > (p.trial.stimulus.timeFpEntered + p.trial.stimulus.fpOffset) % || p.trial.iFrame==p.trial.pldaps.maxFrames) % still need max frames thing here? because we are not ending the trial successfully here anymore
                
                % Turn off fixation point
                p.trial.stimulus.showFixationPoint = 0;
                
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.FIXATION,p.trial.pldaps.iTrial)
                end
                p.trial.ttime      = GetSecs - p.trial.trstart;
                p.trial.stimulus.timeFpOff  = p.trial.ttime;
                p.trial.stimulus.frameFpOff = p.trial.iFrame;
                
                p.trial.state = p.trial.stimulus.states.CHOOSETARG;  
                p.trial.stimulus.showTargets = 1; % probably not needed now
                
            elseif ~fixating && p.trial.ttime < p.trial.stimulus.timeFpEntered + p.trial.stimulus.fpOffset                
                if p.trial.datapixx.use
                    pds.datapixx.flipBit(p.trial.event.BREAKFIX,p.trial.pldaps.iTrial)
                end
                p.trial.stimulus.timeBreakFix = GetSecs - p.trial.trstart;
                p.trial.state = p.trial.stimulus.states.BREAKFIX;
           
            end
        end
                
end    


% Check target fixation
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
            % turn off stimulus
            p.trial.stimulus.colorFixDot        = p.trial.display.clut.bg;            % fixation point 1 color
            p.trial.stimulus.colorFixWindow     = p.trial.display.clut.bg;           % fixation window color
            p.trial.pldaps.goodtrial = 0;
            p.trial.targOn = 2;
            if p.trial.sound.use && ~isnan( p.trial.stimulus.timeFpEntered) 
                PsychPortAudio('Start', p.trial.sound.breakfix, 1, [], [], GetSecs + .1);
            end
            
%             if p.trial.ttime > p.trial.stimulus.timeBreakFix + p.trial.stimulus.breakFixPenalty
                p.trial.flagNextTrial = true;
%             end
        end
        
        
end
    
    