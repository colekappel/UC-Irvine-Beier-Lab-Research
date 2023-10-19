clear

mouse5507=load('5507cellInBin.mat');
mouse5508=load('corrected5508_spatialMapStuff.mat');

mouse5507LR=load('5507LRCorrect.mat');
mouse5508LR=load('5508LRCor.mat');

% now want to divide each column by the time the mouse spent in each bin to
% get the event rate for each cell in each bin 

%Also want to convert the final matrices to tables and add the
%cell names on the left row bc cells arent in order (they are in increasing
%order but some cells are skipped bc they had no events

FR5507=14.255;
FR5508=14.216;

% event per bin matrices now including left and right directions
cellEventPerBinMatx5507 = mouse5507.spikesForCellInBin;
cellEventPerBinMatx5508 = mouse5508.spikesForCellInBin;

cellEventPerBinMatx5507L=mouse5507LR.spikesForCellInBinL;
cellEventPerBinMatx5507R=mouse5507LR.spikesForCellInBinR;

cellEventPerBinMatx5508L=mouse5508LR.spikesForCellInBinL;
cellEventPerBinMatx5508R=mouse5508LR.spikesForCellInBinR;

% mouse frames in bin matrix now need to include left and right directions
mouseFramesInBin5507 = mouse5507.timeInBinX; %these are in frames have to div by FR to get in secs
mouseFramesInBin5508 = mouse5508.timeInBinX; %these are in frames have to div by FR to get in secs

mouseFramesInBin5507L = mouse5507LR.timeInBinXL; %these are in frames have to div by FR to get in secs
mouseFramesInBin5507R = mouse5507LR.timeInBinXR; %these are in frames have to div by FR to get in secs

mouseFramesInBin5508L = mouse5508LR.timeInBinXL; %these are in frames have to div by FR to get in secs
mouseFramesInBin5508R = mouse5508LR.timeInBinXR; %these are in frames have to div by FR to get in secs

% mouse seconds in bin now need to include left and right directions
mouseSecsInBin5507 = mouseFramesInBin5507/FR5507;
mouseSecsInBin5508 = mouseFramesInBin5508/FR5508;

mouseSecsInBin5507L = mouseFramesInBin5507L/FR5507;
mouseSecsInBin5507R = mouseFramesInBin5507R/FR5507;

mouseSecsInBin5508L = mouseFramesInBin5508L/FR5508;
mouseSecsInBin5508R = mouseFramesInBin5508R/FR5508;

% event rate per bin matrix now need to include left and right directions
%eventRatePerBinMatx5507=cellEventPerBinMatx5507 ./ mouseSecsInBin5507;
%eventRatePerBinMatx5508=cellEventPerBinMatx5508 ./ mouseSecsInBin5508;

eventRatePerBinMatx5507L=cellEventPerBinMatx5507L ./ mouseSecsInBin5507L;
eventRatePerBinMatx5507R=cellEventPerBinMatx5507R ./ mouseSecsInBin5507R;

eventRatePerBinMatx5507=[eventRatePerBinMatx5507L,flip(eventRatePerBinMatx5507R,2)];

eventRatePerBinMatx5508L=cellEventPerBinMatx5508L ./ mouseSecsInBin5508L;
eventRatePerBinMatx5508R=cellEventPerBinMatx5508R ./ mouseSecsInBin5508R;

eventRatePerBinMatx5508=[eventRatePerBinMatx5508L,flip(eventRatePerBinMatx5508R,2)];

% probability from formula now need to expand to left and right direc.
%Pi5507=[];
%Pi5508=[];

Pi5507L=[];
Pi5507R=[];

Pi5508L=[];
Pi5508R=[];

%code to get probability mouse in ith bin
for i =1:10
    %Pi5507(i)=mouseSecsInBin5507(i)/(sum(mouseSecsInBin5507));
    %Pi5508(i) =mouseSecsInBin5508(i)/(sum(mouseSecsInBin5508));

    Pi5507L(i)=mouseSecsInBin5507L(i)/(sum(mouseSecsInBin5507L));
    Pi5507R(i)=mouseSecsInBin5507R(i)/(sum(mouseSecsInBin5507R));

    Pi5508L(i) =mouseSecsInBin5508L(i)/(sum(mouseSecsInBin5508L));
    Pi5508R(i) =mouseSecsInBin5508R(i)/(sum(mouseSecsInBin5508R));
end
Pi5507 =[Pi5507L,flip(Pi5507R,2)];
Pi5508 =[Pi5508L,flip(Pi5508R,2)];

% info score now need to do for left and right directions
infoScore5507=infoScore(eventRatePerBinMatx5507, Pi5507);
infoScore5508=infoScore(eventRatePerBinMatx5508, Pi5508);

%infoScore5507L=infoScore(eventRatePerBinMatx5507L, Pi5507L);
%infoScore5507R=infoScore(eventRatePerBinMatx5507R, Pi5507R);

%infoScore5508L=infoScore(eventRatePerBinMatx5508L, Pi5508L);
%infoScore5508R=infoScore(eventRatePerBinMatx5508R, Pi5508R);


%% vvv Code to shuffle spike matrix 100 times vvv now need to do for LR direc.

%vvv start w 5507 vvv this code gets the rando spike matx
infoScoreMatx5507=[];
for i = 1:100
    randoSpikeMatx5507 = randoSpikeMatrix(cellEventPerBinMatx5507);
    randoEventRatePerBinMatx5507=randoSpikeMatx5507 ./ mouseSecsInBin5507; 
    randoInfoScore5507 = infoScore(randoEventRatePerBinMatx5507, Pi5507);
    infoScoreMatx5507(:,i)=transpose(randoInfoScore5507);
end

infoScoreMatx5508=[];
for i = 1:100
    randoSpikeMatx5508 = randoSpikeMatrix(cellEventPerBinMatx5508);
    randoEventRatePerBinMatx5508=randoSpikeMatx5508 ./ mouseSecsInBin5508; 
    randoInfoScore5508 = infoScore(randoEventRatePerBinMatx5508, Pi5508);
    infoScoreMatx5508(:,i)=transpose(randoInfoScore5508);
end

%now loop through all cells in the original info score matrix 
% and mark it as a place cell if its value is higher 
% than 95th percentile of its corresponding row 
% now need to do for LR direc.s + want to plot as heat maps
placeCells5507="";
cellNames5507 = table2array(mouse5507.cellNames);
b=1;
for i = 1:length(infoScore5507)
    if infoScore5507(i)>prctile(infoScoreMatx5507(i,:),95)
        placeCells5507(b)=string(cellNames5507(i));
        placeCells5507Arr(b,:)=eventRatePerBinMatx5507(i,:);
        infoScoresplaceCells5507(b)=infoScore5507(i);
        b=b+1;
    end
end

placeCells5508="";
cellNames5508 = table2array(mouse5508.cellNames);
b=1;
for i = 1:length(infoScore5508)
    if infoScore5508(i)>prctile(infoScoreMatx5508(i,:),95)
        placeCells5508(b)=string(cellNames5508(i));
        placeCells5508Arr(b,:)=eventRatePerBinMatx5508(i,:);
        infoScoresplaceCells5508(b)=infoScore5508(i);
        b=b+1;
    end
end

%now code to visualize diagonal line: start w 5507
orderMatx5507=[];

for i =1:length(infoScoresplaceCells5507)
    [m,j]=max(placeCells5507Arr(i,:));
    orderMatx5507(i)=j;
end

[B,I]=sort(orderMatx5507);
sortedPlaceCells5507=[];
for i = 1: length(orderMatx5507)
    sortedPlaceCells5507(i,:)=placeCells5507Arr(I(i),:);
end
heatmap(sortedPlaceCells5507);
title("5507 Place Cells in Order of Bin")
xlabel("Bin (1-10 left moving bin, 11-20 right moving bin)")
ylabel("Cell")
figure;

orderMatx5508=[];

for i =1:length(infoScoresplaceCells5508)
    [m,j]=max(placeCells5508Arr(i,:));
    orderMatx5508(i)=j;
end

[B,I]=sort(orderMatx5508);
sortedPlaceCells5508=[];
for i = 1: length(orderMatx5508)
    sortedPlaceCells5508(i,:)=placeCells5508Arr(I(i),:);
end
heatmap(sortedPlaceCells5508);
title("5508 Place Cells in Order of Bin")
xlabel("Bin (1-10 left moving bin, 11-20 right moving bin)") %11 corresponds to bin 10, 12 to bin 9 ... 20 to bin 1
ylabel("Cell")

save('Place Cell Analysis.mat')