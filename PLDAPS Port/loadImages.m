function dv = loadImages(dv)
% dv.filePaths is now set in liprein (the condition file)

switch dv.filePaths
    case 'testingStim'
        % Sets Paths
        dv.fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\novel_objects\png_files\';
        dv.fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\natural\';
        
        % File formats
        dv.fileInfo.objectFormat = 'png';
        dv.fileInfo.sceneFormat = 'bmp';
        
    case 'conservativeStim'
        % Sets Paths
        dv.fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\shapes\bow\';
        dv.fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\natural\';
        
        % File formats
        dv.fileInfo.objectFormat = 'tif';
        dv.fileInfo.sceneFormat = 'bmp';
        
    case 'naturalisticStim'
        % Sets Paths
        dv.fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\shapes\wob\';
        dv.fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\baboon_habitat\resized\newAR\';
        
        % File formats
        dv.fileInfo.objectFormat = 'tif';
        dv.fileInfo.sceneFormat = 'jpg';
        
    case 'realDeal'
        % Sets Paths
        dv.fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\shapes\bow\';
        dv.fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\lipreinst_stimuli\scenemodel\natural40\';
        
        % File formats
        dv.fileInfo.objectFormat = 'tif';
        dv.fileInfo.sceneFormat = 'jpg';
        
    case 'shopRig'
        % Sets Paths
        dv.fileInfo.objectPath = '/Users/huklab/ELH/Stimuli/shapes/bow/';  % different slash on Macs
        dv.fileInfo.scenePath = '/Users/huklab/ELH/Stimuli/natural40/';
        
        % File formats
        dv.fileInfo.objectFormat = 'tif';
        dv.fileInfo.sceneFormat = 'jpg';
        
end
% Get Images from directories
objectImages = dir([dv.fileInfo.objectPath '*.' dv.fileInfo.objectFormat]); 
sceneImages = dir([dv.fileInfo.scenePath '*.' dv.fileInfo.sceneFormat]);

% Randomize Images and pair them
% objectImages_mask = Shuffle(1:numel(objectImages));
% sceneImages_mask = Shuffle(1:numel(sceneImages));
% objectImages_mask = Shuffle(repmat(objectImages_mask, 1, length(sceneImages_mask) / length(objectImages_mask)));

% Testing randperm instead because Shuffle depends on psychtoolbox (?)
% sceneImages_mask = randperm(length(sceneImages));
% objectImages_mask = datasample(repmat(randperm(length(objectImages)),1,(length(sceneImages) / length(objectImages))),length(sceneImages),'Replace',false);

% Eliminating use of datasample because this might be having an aversive effect on randomization
sceneImages_mask = randperm(length(sceneImages));
objectImages_mask = repmat(randperm(length(objectImages)),1,(length(sceneImages) / length(objectImages)));
objectImages_mask = objectImages_mask(randperm(length(objectImages_mask)));

% Saving order
dv.pairOrder = cell(numel(objectImages_mask),4);
for i = 1:length(dv.pairOrder)
    dv.pairOrder{i,1} = objectImages(objectImages_mask(i)).name;
    dv.pairOrder{i,2} = sceneImages(sceneImages_mask(i)).name;
    dv.pairOrder{i,3} = objectImages_mask(i);
    dv.pairOrder{i,4} = sceneImages_mask(i);
end

% Proliferate pairings for many trials - 40 is the current size of a
% stimulus set: 8 shapes, each paired with 5 scenes
% for i = 1:dv.finish / size(dv.pairOrder,1) 
%     dv.pairOrder = [dv.pairOrder;
%     datasample(dv.pairOrder,size(dv.pairOrder,1),1,'replace',false)]; % should we overwrite dv.pairOrder, this was the problem!!!!!!
% end

pairOrderTemp = cell(numel(dv.finish),4);
pairOrderTemp = datasample(dv.pairOrder,size(dv.pairOrder,1),1,'replace',false);
for i = 1:dv.finish / size(dv.pairOrder,1)
    pairOrderRep = datasample(dv.pairOrder,size(dv.pairOrder,1),1,'replace',false);
    pairOrderTemp = vertcat(pairOrderTemp, pairOrderRep);
end

dv.pairOrder = pairOrderTemp;

% Attach Session ID
sessionID = repmat(cellstr(dv.pref.sfile), length(dv.pairOrder), 1); % is dv.finish or the +40 going to be a problem here?
dv.pairOrder = horzcat(dv.pairOrder, sessionID);

end

