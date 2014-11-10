function [Dx, Dy] = calcVisAngDS(alpha, Z, scrX, scrY, screenXpixels, screenYpixels)
% Visual Angle calculations
% [Dx, Dy] = calcVisAng(alpha, Z, scrX, scrY)  % physical screen size in
% cm, Z in cm (dist. to screen), alpha in degrees

% alpha = 6; % Visual angle in degrees (it will be turned into radians)
% Z = 57; % Distance from eyes to the screen (in cm?)
% scrX = 16 * 2.54; % real screen size in inches, convert to cm
% scrY = 9 * 2.54;
Dx = Z * screenXpixels / (scrX) * tan(alpha*(pi/180)); % Get pixel dimensions from psychtoolbox call, also uses pixels per cm here, necessary?
Dy = Z * screenYpixels / (scrY) * tan(alpha*(pi/180));
end


% notes for calculating with dependencies between x and y
% distance in x and y in pixels first?? need those measurements first
% Dx = sqrt(Z^2 + y^2) * screenXpixels / (scrX * 2.54) * tan(alpha*(pi/180));
% Dy = sqrt(Z^2 + x^2) * screenYpixels / (scrY * 2.54) * tan(alpha*(pi/180));