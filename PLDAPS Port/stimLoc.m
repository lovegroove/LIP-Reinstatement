function dv = stimLoc(dv)
% calculates the location of the stimulus based on angle from the center and distance
% needs dv.trial.theta input for desired angle
% also relies upon calcVisAngDS function and screen parameters

r = sqrt((dv.pa.Dx)^2 + (dv.pa.Dy)^2); % degrees of visual angle (which r will be the x and y dimensions)
dv.trial.object1loc = [dv.pa.xCenter + r * cos(dv.trial.theta) dv.pa.yCenter + r * sin(dv.trial.theta)];
dv.trial.object2loc = [dv.pa.xCenter + r * cos(dv.trial.theta+pi) dv.pa.yCenter + r * sin(dv.trial.theta+pi)];

end