function stimMatrix=buildStimMatrix(stimSequence, stimIntervals, time)

stimMatrix = false(length(unique(stimSequence)), length(time));

for iStim = 1:length(stimSequence)
    ind = (time >= stimIntervals(iStim,1) & time <= stimIntervals(iStim,2));
    stimMatrix(stimSequence(iStim), ind) = true;
end



