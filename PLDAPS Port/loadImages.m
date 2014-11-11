function dv = loadImages(dv)
% [sceneOrder, objectOrder, objectPath, scenePath, objectFormat, sceneFormat, objectImages_mask, sceneImages_mask, objectImages, sceneImages] = loadImagesFunc

% dv.filePaths is now set in liprein (the condition file)
%dv.filePaths = 'shopRig'; % testingStim, conservativeStim, or naturalisticStim, for shapes - bow or wob

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
objectImages_mask = Shuffle(1:numel(objectImages));
sceneImages_mask = Shuffle(1:numel(sceneImages));
objectImages_mask = Shuffle(repmat(objectImages_mask, 1, length(sceneImages_mask) / length(objectImages_mask)));

% Saving order
dv.pairOrder = cell(numel(objectImages_mask),4);
for i = 1:length(dv.pairOrder)
    dv.pairOrder{i,1} = objectImages(objectImages_mask(i)).name;
    dv.pairOrder{i,2} = sceneImages(sceneImages_mask(i)).name;
    dv.pairOrder{i,3} = objectImages_mask(i);
    dv.pairOrder{i,4} = sceneImages_mask(i);
end

end

