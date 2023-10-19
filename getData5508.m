clear

spikes5508 = readtable('Data to Analyze/5508.spikes.LT-imaging-T1.csv');

cellNames=unique(spikes5508(:,2)); %array of unique cell names in spike data

% use ismember instead of == to compare 2 table elements (same as ==)

spikeDataLength=height(spikes5508);
%try for 1 cell first:
i=1;
avAmplitudeForEachCell =[];
amplitudeArray=[];
spikesPerMouse=[];
numCells=height(cellNames);
for p = 1:numCells
     z=1;
     a=0;
     while ismember(cellNames(p,1),spikes5508(i,2))
        %newTable(i,:)=spikes5507(i,:);
        amplitudeArray(z)=spikes5508{i,3};
        i=i+1;
        z=z+1;
        a=a+1;
        if i>spikeDataLength
            break
        end
     end
     spikesPerMouse(p)=a;
     avAmplitudeForEachCell(p)=mean(amplitudeArray);
     amplitudeArray=[];
end

avAmplitudeForEachCell=transpose(avAmplitudeForEachCell);

cellWAvAmp= [cellNames array2table(avAmplitudeForEachCell)];

totalTimeSpikesRecorded = 995.86-3.869; %specific to 5508 (accurate but not sup precise) in seconds

eventRate = spikesPerMouse/totalTimeSpikesRecorded; %in hertz

eventRate=transpose(eventRate);

cellWAvAmpAndEventRate=[cellWAvAmp array2table(eventRate)];

boxplot(avAmplitudeForEachCell)
title('Average Amplitude for Each Cell in 5508 5xFAD-')