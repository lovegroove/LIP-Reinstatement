function [pairOrder, fileInfo] = loadImagePairsFunc
% [sceneOrder, objectOrder, objectPath, scenePath, objectFormat, sceneFormat, objectImages_mask, sceneImages_mask, objectImages, sceneImages] = loadImagesFunc
filePaths = 'realDeal'; % testingStim, conservativeStim, or naturalisticStim, for shapes - bow or wob

switch filePaths
    case 'testingStim'
        % Sets Paths
        fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\novel_objects\png_files\';
        fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\natural\';
        
        % File formats
        fileInfo.objectFormat = 'png';
        fileInfo.sceneFormat = 'bmp';
        
    case 'conservativeStim'
        % Sets Paths
        fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\shapes\bow\';
        fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\natural\';
        
        % File formats
        fileInfo.objectFormat = 'tif';
        fileInfo.sceneFormat = 'bmp';
        
    case 'naturalisticStim'
        % Sets Paths
        fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\shapes\wob\';
        fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\baboon_habitat\resized\newAR\';
        
        % File formats
        fileInfo.objectFormat = 'tif';
        fileInfo.sceneFormat = 'jpg';
        
    case 'realDeal'
        % Sets Paths
        fileInfo.objectPath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\shapes\bow\';
        fileInfo.scenePath = 'C:\Users\MrLovegroove\Documents\MATLAB\stimuli\lipreinst_stimuli\scenemodel\natural40\';
        
        % File formats
        fileInfo.objectFormat = 'tif';
        fileInfo.sceneFormat = 'jpg';
end
% Get Images from directories
objectImages = dir([fileInfo.objectPath '*.' fileInfo.objectFormat]); 
sceneImages = dir([fileInfo.scenePath '*.' fileInfo.sceneFormat]);

% Randomize Images and pair them
objectImages_mask = Shuffle(1:numel(objectImages));
sceneImages_mask = Shuffle(1:numel(sceneImages));
objectImages_mask = Shuffle(repmat(objectImages_mask, 1, length(sceneImages_mask) / length(objectImages_mask)));

% Saving order
pairOrder = cell(numel(objectImages_mask),4);
for i = 1:length(pairOrder)
    pairOrder{i,1} = objectImages(objectImages_mask(i)).name;
    pairOrder{i,2} = sceneImages(sceneImages_mask(i)).name;
    pairOrder{i,3} = objectImages_mask(i);
    pairOrder{i,4} = sceneImages_mask(i);
end

end