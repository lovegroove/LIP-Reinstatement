% Save Pairings in mat file
saveFile = [dv.subj date];
save(saveFile,'-struct','dv','pairOrder');