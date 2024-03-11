function[] = AnalysisFunction(spikes,ezTrack,FR,inscopixStartTimeInSeconds,ezTrackStart,...
    ezTrackEnd,mouseNum,PLUSORMINUS,NUMBEROFBINS,cohortID)
% use 1 for plus, 0 for minus ie:
% vvv these need to be edited for each mouse vvv 
%{
spikes = readtable('Data to Analyze/5508.spikes.LT-imaging-T1.csv');
ezTrack = readtable('Data to Analyze/WIN_20230601_5508_LT_Imaging_T1_LocationOutput.csv');
FR=14.216;
inscopixStartTimeInSeconds = 55/FR; %frames over frame rate
ezTrackStart=270/FR; %frames over frame rate
ezTrackEnd=13804/FR; %frames over frame rate
mouseNum = "5508";
PLUSORMINUS=0; % 1 for plus 0 for minus
%}
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



ezTrack{:,12}=ezTrack{:,7}/FR; % add time in seconds to table in 12th column
ezTrack.Properties.VariableNames([12])="Time in Seconds";

%cnmfe5508 = readmatrix('Data to Analyze/5508.cnmfe.LT-imaging-T1.csv');

% get the time of the spikes 
timeOfSpikes = table2array(spikes(1:end,1:1));

%vvvv code to get the ez track and inscopix times aligned vvvv
inscopixStartTimesToDelete = ezTrackStart-inscopixStartTimeInSeconds;
inscopixEndTimesToDelete=ezTrackEnd-inscopixStartTimeInSeconds;
indices = find(timeOfSpikes<inscopixStartTimesToDelete);
indicesEndOfVid=find(timeOfSpikes>inscopixEndTimesToDelete);
indices = [indices;indicesEndOfVid];

spikesCorrected = spikes;
spikesCorrected(indices,:)=[];

%shift all values to be aligned with ez track
spikesCorrected(:,1) = spikesCorrected(:,1)-inscopixStartTimesToDelete;

min(spikesCorrected(:,1))

spikeTimeCorrected = table2array(spikesCorrected(:,1));
lenSpikeTime = height(spikeTimeCorrected);
lenEzTrack = height(ezTrack);
xValMatched = zeros(lenSpikeTime,1);

allXVals = table2array(ezTrack(:,8));
allTimesEzTrack = table2array(ezTrack(:,12));

for i = 1:lenSpikeTime

    for p = 1:lenEzTrack
     % want to say if abs val of diff between times less than 0.05 then take
        % that x val and break the loop too

        if abs(spikeTimeCorrected(i)-allTimesEzTrack(p))<0.05
            %then add the x val to the xVal array
            xValMatched(i)=allXVals(p);
        end

    end

end

xValMatchedTable = array2table(xValMatched,'VariableNames',{'X Value'});

spikesWXVal = [spikesCorrected xValMatchedTable];

%add in line to delete values greater than the max time in allTimesEzTrack
moreIndices=find(spikesWXVal{:,1}>max(allTimesEzTrack));
spikesWXVal(moreIndices,:)=[];

%% vvvvvvvvvvvvv Code to delete end bins vvvvvvvvvvvvvvvvv
%NUMBEROFBINS=10; %change to desired # of bins
BINSIZE = max(spikesWXVal{:,4})/NUMBEROFBINS; %10 needs to be changed to number of bins
RowsToDelete=find(spikesWXVal{:,4}<BINSIZE);
spikesWXVal(RowsToDelete,:)=[];
MoreRowsToDelete=find(spikesWXVal{:,4}>(BINSIZE*(NUMBEROFBINS-1)));
spikesWXVal(MoreRowsToDelete,:)=[];

xvalRTD=find(allXVals<BINSIZE);
allXVals(xvalRTD)=[];
allTimesEzTrack(xvalRTD)=[];
morexvalRTD=find(allXVals>(BINSIZE*(NUMBEROFBINS-1)));
allXVals(morexvalRTD)=[];
allTimesEzTrack(morexvalRTD)=[];
%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%writetable(spikesWXVal5508,'spikesWXVal5508.csv');

%Now want to make a dv/dt array so we can label mouse as moving left or
%right

%this for loop will label xvals starting from 2 to (end - 1) th point
%if its not moving we say its moving left arbitrarily and also at end
%points say its moving left arbitrarily

leftOrRight="";
for i = 1:length(allXVals)
    if i ==1 || i==length(allXVals)
        leftOrRight(i)='L';
    elseif allXVals(i+1)-allXVals(i-1)>0
        leftOrRight(i)='R';
    elseif allXVals(i+1)-allXVals(i-1)<=0
        leftOrRight(i)='L';
    end
end

leftOrRight=transpose(leftOrRight);
allXValsWDirection = [allXVals,leftOrRight];

%now want to attach left or right direction to spike data based on the
%spike point in time
leftOrRight_spikes=strings([1,height(spikesWXVal)]); 
%^^ that has L or R for spikeWXVal5508
for i = 1:height(spikesWXVal)
    %disp(i);
    %disp("Line 95");
    [row,col]=find(allTimesEzTrack>spikesWXVal{i,1}-.5);
    for j = row(1):row(end)
        if (spikesWXVal{i,1}-allTimesEzTrack(j))<=0.08
            leftOrRight_spikes(i)=leftOrRight(j);
            break
        end
    end
end
cellNames=unique(spikesWXVal(:,2)); %array of unique cell names in spike data
leftOrRight_spikes=transpose(leftOrRight_spikes);
spikesWXVal.LeftOrRight = leftOrRight_spikes;

[Pi,mouseSecsInBinX,eventPerBinMatx,eventRatePerBinMatx]...
    = getEventRatePerBinMatx(NUMBEROFBINS,BINSIZE,allXValsWDirection,cellNames,spikesWXVal,leftOrRight_spikes,FR);

theInfoScore=infoScore(eventRatePerBinMatx, Pi);

%vvv get sparsity and spatial coherence vvv

sparsity=zeros(1,height(eventRatePerBinMatx));
for i = 1: height(eventRatePerBinMatx)
    top = 0;
    bottom=0;

    for j = 1:width(eventRatePerBinMatx)
        top=top+(Pi(j)*eventRatePerBinMatx(i,j))^2;
        bottom=bottom+Pi(j)*(eventRatePerBinMatx(i,j)^2);
    end
    if bottom ==0
        sparsity(i)=0;
    else
    sparsity(i)=top/bottom;
    end
end

%calculate spatial coherence, start w 5507

surroundingAv=zeros(height(eventRatePerBinMatx),width(eventRatePerBinMatx));
%want to calculate the average of the surrounding bins
for i = 1:height(eventRatePerBinMatx)
    for j=1:width(eventRatePerBinMatx)
        if j==1
            surroundingAv(i,j)=eventRatePerBinMatx(i,j+1);
        elseif j==width(eventRatePerBinMatx)
            surroundingAv(i,j)=eventRatePerBinMatx(i,j-1);
        else
            surroundingAv(i,j)=(eventRatePerBinMatx(i,j-1)+eventRatePerBinMatx(i,j+1))/2;
        end
    end
end

spatialCoherence = zeros(1,height(eventRatePerBinMatx));
for i = 1:height(eventRatePerBinMatx)
   spatialCoherence(i)=corr(transpose(eventRatePerBinMatx(i,1:width(eventRatePerBinMatx))),transpose(surroundingAv(i,1:width(eventRatePerBinMatx))),'type','Pearson');
      if isnan(spatialCoherence(i))
        spatialCoherence(i)=0;
    end
end
spatialCoherence=spatialCoherence.^2;

%code to get place cells
infoScoreMatx=[];
for i = 1:100
    randoSpikeMatx = randoSpikeMatrix(eventPerBinMatx); 
    randoEventRatePerBinMatx=randoSpikeMatx ./ mouseSecsInBinX; 
    randoInfoScore = infoScore(randoEventRatePerBinMatx, Pi);
    infoScoreMatx(:,i)=transpose(randoInfoScore);
end

%now loop through all cells in the original info score matrix 
% and mark it as a place cell if its value is higher 
% than 95th percentile of its corresponding row 
% now need to do for LR direc.s + want to plot as heat maps
placeCells="";
cellNames=table2array(unique(spikesWXVal(:,2)));
b=1;
for i = 1:length(theInfoScore)
    if theInfoScore(i)>prctile(infoScoreMatx(i,:),95)
        placeCells(b)=string(cellNames(i));
        placeCellsArr(b,:)=eventRatePerBinMatx(i,:);
        infoScoresplaceCells(b)=theInfoScore(i);
        b=b+1;
    end
end

%now code to visualize diagonal line: start w 5507
orderMatx=[];

for i =1:length(infoScoresplaceCells)
    [m,j]=max(placeCellsArr(i,:));
    orderMatx(i)=j;
end

[B,I]=sort(orderMatx);
sortedPlaceCells=[];
for i = 1: length(orderMatx)
    sortedPlaceCells(i,:)=placeCellsArr(I(i),:);
end

%% code to get av amp per cell and event rate: copy pasted from other program
spikeDataLength=height(spikesWXVal);
totalTimeSpikesRecorded=max(spikesWXVal{:,1})-min(spikesWXVal{:,1});
%try for 1 cell first:
i=1;
avAmplitudeForEachCell =[];
amplitudeArray=[];
spikesPerCell=[];
cellNamesNew=unique(spikesWXVal(:,2));
numCells=height(cellNamesNew);
for p = 1:numCells
     z=1;
     a=0;
     sumAmplitude=0;
     while ismember(cellNamesNew(p,1),spikesWXVal(i,2))
        %newTable(i,:)=spikes5507(i,:);
        amplitudeArray(z)=spikesWXVal{i,4};
        i=i+1;
        z=z+1;
        a=a+1;
        if i>spikeDataLength
            break
        end
     end
     spikesPerCell(p)=a;
     avAmplitudeForEachCell(p)=mean(amplitudeArray);
     amplitudeArray=[];
end

eventRatePerCell = spikesPerCell/totalTimeSpikesRecorded; %in hertz

%% want to make tables for linear effects models:
% For Amplitude:
avAmplitudeForEachCell=transpose(avAmplitudeForEachCell);

if PLUSORMINUS==1
    FadPlusorMinus = ones(length(avAmplitudeForEachCell),1); %1 for plus, 0 for minus
elseif PLUSORMINUS==0
    FadPlusorMinus = zeros(length(avAmplitudeForEachCell),1); %1 for plus, 0 for minus
end
MouseID=strings(length(avAmplitudeForEachCell),1);
for i = 1:length(avAmplitudeForEachCell)
    MouseID(i)=mouseNum;
end
LMETableAmplitude=table(avAmplitudeForEachCell,FadPlusorMinus,MouseID);

%for event rate:
eventRatePerCell = transpose(eventRatePerCell);
LMETableEventRate=table(eventRatePerCell,FadPlusorMinus,MouseID);

%for info score
InformationScore=transpose(theInfoScore);
LMETableInfoScore=table(InformationScore,FadPlusorMinus,MouseID);

%for sparsity
sparsity=transpose(sparsity);
LMETableSparsity=table(sparsity,FadPlusorMinus,MouseID);

%for spatial coherence
spatialCoherence=transpose(spatialCoherence);
LMETableSpatC=table(spatialCoherence,FadPlusorMinus,MouseID);

%% code to get stability 
num = floor(length(eventRatePerCell)/2);
pearsonCorrER = corr(eventRatePerCell(1:num),eventRatePerCell((num+1):(num*2)),'type','Pearson');
stability = atan(pearsonCorrER);

%% code to get behavior data

% create derivative array
derivative=zeros(1,length(allXVals)); % x coord / time
for j = 2:(length(allXVals)-1)
    derivative(j)=(allXVals(j+1)-allXVals(j-1))/(allTimesEzTrack(j+1)-allTimesEzTrack(j-1));
end
turnLocation=[]; %will be where derivative is zero
turnTime=[];
threshold =1; %for finding where turns are 

crossingCords=[];
timePerCrossing=[];
g=1;
iCord=[];
p=1;
i=1;
timePerCrossing_Time=[];
timePerCrossing_X=[];
qqq=1;
zzz=1;
while i< length(allXVals)
        iCord(p)=i;
        p=p+1;

    if allXVals(i)>=700 && allXVals(i)<=950
        index=i;
        crossingCords(g)=allXVals(i); %i is the index
        timeStampForCrossingCords(g)=allTimesEzTrack(i);
  
        while allXVals(i) <=1340 && allXVals(i) >=200
            if i>=length(derivative)
                break
            end
            i=i+1;
            
            if derivative(i)<=threshold && derivative(i)>=-threshold
                turnLocation(zzz)=allXVals(i);
                turnTime(zzz)=allTimesEzTrack(i);
                zzz=zzz+1;
            end
        end
        endIndex=i;

        while allXVals(index) >=200 && allXVals(index) <=1340
            if index==1
                break
            end
            index=index-1;
            if derivative(index)<=threshold && derivative(index)>=-threshold
                turnLocation(zzz)=allXVals(index);
                turnTime(zzz)=allTimesEzTrack(index);
                zzz=zzz+1;
            end
        end

        startIndex=index;
        timePerCrossing(g)=abs(allTimesEzTrack(endIndex)-allTimesEzTrack(startIndex));
        distanceTraveledPerCrossing(g)=abs(allXVals(endIndex)-allXVals(startIndex));
        speedPerCrossing(g)= distanceTraveledPerCrossing(g)/timePerCrossing(g);
        g=g+1;
        if index>=i
            i=index;
        end
        % vvvvv just for visual 2 vvvvvvvv
        timePerCrossing_X(qqq)=allXVals(index);
        timePerCrossing_Time(qqq)=allTimesEzTrack(index);
        qqq=qqq+1;
        timePerCrossing_X(qqq)=allXVals(i);
        timePerCrossing_Time(qqq)=allTimesEzTrack(i);
        qqq=qqq+1;
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    else 
        i=i+1;
    end
end
averageSpeedPerCrossing=mean(speedPerCrossing);
numberOfCrossings= length(crossingCords);

[Placeholder,Ind] =sort(turnTime);
turnTime=turnTime(Ind);
turnLocation = turnLocation(Ind);

deletions = find(diff(turnTime)<2);
turnTime(deletions)=[];
turnLocation(deletions)=[];

numTurnsNotAtEnds=length(turnLocation);

%plot data to ensure accuracy of the method:
%{
plot(allTimesEzTrack,allXVals)
hold on 
plot(timeStampForCrossingCords,crossingCords,'.')
hold on
plot(timePerCrossing_Time,timePerCrossing_X,'o','Color','r')
hold on
plot(turnTime,turnLocation,'o','Color','magenta')
%}

%% calculate stability: first will look at 1st half of time vs 2nd half of time
% first split allXVals into 1st and 2nd half
[RowsAllX1,ColsAllx1]=find(allTimesEzTrack<=(max(allTimesEzTrack)/2));
[RowsAllX2,ColsAllx2]=find(allTimesEzTrack>(max(allTimesEzTrack)/2));
allXValsWDirection1 = allXValsWDirection(RowsAllX1,:); %1st half allXValsWDirection
allXValsWDirection2 = allXValsWDirection(RowsAllX2,:); %2nd half allXValsWDirection

%Now split up spikesWXVal and leftOrRight_spikes
[RowsSpikes1,ColsSpikes1]=find(spikesWXVal{:,1}<=(max(allTimesEzTrack)/2));
[RowsSpikes2,ColsSpikes2]=find(spikesWXVal{:,1}>(max(allTimesEzTrack)/2));
spikesWXVal1 =spikesWXVal(RowsSpikes1,:);
spikesWXVal2 =spikesWXVal(RowsSpikes2,:);
leftOrRight_spikes1=leftOrRight_spikes(RowsSpikes1,:);
leftOrRight_spikes2=leftOrRight_spikes(RowsSpikes2,:);

% now we just need to calculate eventRatePerBinMatx for 1st half and 2nd half
% and calculate stability from the formula (atan of pearson correlation)
%1st half:
[Pi1,mouseSecsInBinX1,eventPerBinMatx1,eventRatePerBinMatx1]...
    = getEventRatePerBinMatx(NUMBEROFBINS,BINSIZE,allXValsWDirection1,cellNames,spikesWXVal1,leftOrRight_spikes1,FR);
%2nd half:
[Pi2,mouseSecsInBinX2,eventPerBinMatx2,eventRatePerBinMatx2]...
    = getEventRatePerBinMatx(NUMBEROFBINS,BINSIZE,allXValsWDirection2,cellNames,spikesWXVal2,leftOrRight_spikes2,FR);

stability1v2 = zeros(1,height(eventRatePerBinMatx));
for i = 1:height(eventRatePerBinMatx)
   stability1v2(i)=atan(corr(transpose(eventRatePerBinMatx1(i,1:width(eventRatePerBinMatx))),transpose(eventRatePerBinMatx2(i,1:width(eventRatePerBinMatx))),'type','Pearson')^2);
    if isnan(stability1v2(i))
        stability1v2(i)=0;
    end
end

%% now to calculate stability of odd vs even trials
% start w/ odd trials:
ii = 1;
spikesWXVal11=[];
allXValsWDirection11=[];
leftOrRight_spikes11=[];
while ii < (length(timePerCrossing_Time)-3)
    [rowsSpikes11,colsSpikes11]=...
        find(spikesWXVal{:,1}>=timePerCrossing_Time(ii) & spikesWXVal{:,1}<=timePerCrossing_Time(ii+3));
    [rowsAllX11,colsAllX11]=...
        find(allTimesEzTrack(:,1)>=timePerCrossing_Time(ii) & allTimesEzTrack(:,1)<=timePerCrossing_Time(ii+3));
    spikesWXVal11=[spikesWXVal11;spikesWXVal(rowsSpikes11,:)];
    leftOrRight_spikes11=[leftOrRight_spikes11;leftOrRight_spikes(rowsSpikes11,:)];
    allXValsWDirection11=[allXValsWDirection11;allXValsWDirection(rowsAllX11,:)];

    ii=ii+8;
end

% Now for even trials:
ii = 5;
spikesWXVal22=[];
allXValsWDirection22=[];
leftOrRight_spikes22=[];
while ii < (length(timePerCrossing_Time)-3)
    [rowsSpikes22,colsSpikes22]=...
        find(spikesWXVal{:,1}>=timePerCrossing_Time(ii) & spikesWXVal{:,1}<=timePerCrossing_Time(ii+3));
    [rowsAllX22,colsAllX22]=...
        find(allTimesEzTrack(:,1)>=timePerCrossing_Time(ii) & allTimesEzTrack(:,1)<=timePerCrossing_Time(ii+3));
    spikesWXVal22=[spikesWXVal22;spikesWXVal(rowsSpikes22,:)];
    leftOrRight_spikes22=[leftOrRight_spikes22;leftOrRight_spikes(rowsSpikes22,:)];
    allXValsWDirection22=[allXValsWDirection22;allXValsWDirection(rowsAllX22,:)];

    ii=ii+8;
end

% and now to calculate event rate per bin matx for each

%odd:
[Pi11,mouseSecsInBinX11,eventPerBinMatx11,eventRatePerBinMatx11]...
    = getEventRatePerBinMatx(NUMBEROFBINS,BINSIZE,allXValsWDirection11,cellNames,spikesWXVal11,leftOrRight_spikes11,FR);
%even:
[Pi22,mouseSecsInBinX22,eventPerBinMatx22,eventRatePerBinMatx22]...
    = getEventRatePerBinMatx(NUMBEROFBINS,BINSIZE,allXValsWDirection22,cellNames,spikesWXVal22,leftOrRight_spikes22,FR);

stabilityOddEven = zeros(1,height(eventRatePerBinMatx));
for i = 1:height(eventRatePerBinMatx)
   stabilityOddEven(i)=atan(corr(transpose(eventRatePerBinMatx11(i,1:width(eventRatePerBinMatx))),transpose(eventRatePerBinMatx22(i,1:width(eventRatePerBinMatx))),'type','Pearson')^2);
    if isnan(stabilityOddEven(i))
        stabilityOddEven(i)=0;
    end
end

stabilityAverage=zeros(1,length(stabilityOddEven));
for i=1:length(stabilityOddEven)
    stabilityAverage(i)=mean([stability1v2(i);stabilityOddEven(i)]);
end

%% stability LME Tables
stability1v2=transpose(stability1v2);
LMETableStability1v2=table(stability1v2,FadPlusorMinus,MouseID);

stabilityOddEven=transpose(stabilityOddEven);
LMETableStabilityOE=table(stabilityOddEven,FadPlusorMinus,MouseID);

stabilityAverage=transpose(stabilityAverage);
LMETableStabilityAv=table(stabilityAverage,FadPlusorMinus,MouseID);

%% time per crossing table: different for each mouse
timePerCrossing=transpose(timePerCrossing);
index = length(timePerCrossing)+1; % have to modify arrays so
%FadPlusorMinus(index:end,:)=[]; % they all have same length
%MouseID(index:end,:)=[];
%LMETableTimePerCrossing=table(timePerCrossing,FadPlusorMinus,MouseID);

if length(timePerCrossing)>length(FadPlusorMinus)
    for i = 1:length(timePerCrossing)
        FadPlusorMinus(i)=PLUSORMINUS;
        MouseID(i)=mouseNum;
    end
else 
    FadPlusorMinus(index:end,:)=[]; % they all have same length
    MouseID(index:end,:)=[];
    LMETableTimePerCrossing=table(timePerCrossing,FadPlusorMinus,MouseID);
end

LMETableTimePerCrossing=table(timePerCrossing,FadPlusorMinus,MouseID);

if cohortID==0
    LMETableTimePerCrossing.cohortID = zeros(height(LMETableTimePerCrossing),1);
    LMETableStability1v2.cohortID = zeros(height(LMETableStability1v2),1);
    LMETableStabilityOE.cohortID = zeros(height(LMETableStabilityOE),1);
    LMETableStabilityAv.cohortID = zeros(height(LMETableStabilityAv),1);
    LMETableAmplitude.cohortID = zeros(height(LMETableAmplitude),1);
    LMETableEventRate.cohortID = zeros(height(LMETableEventRate),1);
    LMETableInfoScore.cohortID = zeros(height(LMETableInfoScore),1);
    LMETableSparsity.cohortID = zeros(height(LMETableSparsity),1);
    LMETableSpatC.cohortID = zeros(height(LMETableSpatC),1);

elseif cohortID==1
    LMETableTimePerCrossing.cohortID = ones(height(LMETableTimePerCrossing),1);
    LMETableStability1v2.cohortID = ones(height(LMETableStability1v2),1);
    LMETableStabilityOE.cohortID = ones(height(LMETableStabilityOE),1);
    LMETableStabilityAv.cohortID = ones(height(LMETableStabilityAv),1);
    LMETableAmplitude.cohortID = ones(height(LMETableAmplitude),1);
    LMETableEventRate.cohortID = ones(height(LMETableEventRate),1);
    LMETableInfoScore.cohortID = ones(height(LMETableInfoScore),1);
    LMETableSparsity.cohortID = ones(height(LMETableSparsity),1);
    LMETableSpatC.cohortID = ones(height(LMETableSpatC),1);
end

filename = strcat("Mouse Data/Data_",mouseNum,"EndBinsRemoved");
disp(filename)
save(filename,"stabilityAverage","stability1v2","stabilityOddEven",...
    "theInfoScore","sparsity","spatialCoherence","avAmplitudeForEachCell",...
    "eventRatePerCell","LMETableTimePerCrossing","timePerCrossing",...
    "LMETableStability1v2","LMETableStabilityOE","LMETableStabilityAv",...
    "LMETableAmplitude","LMETableEventRate","LMETableInfoScore",...
    "LMETableSparsity","LMETableSpatC","placeCellsArr","infoScoresplaceCells",...
    "spikesWXVal","placeCells")
end