clear

%Note: timeInBinX array is in frames, have to div by FR to get in seconds
%% first part of code
spikes5507 = readtable('Data to Analyze/5507.spikes.LT-imaging-T1.csv');
ezTrack5507 = readtable('Data to Analyze/WIN_20230601_5507_LT_Imaging_T1_LocationOutput.csv');
%cnmfe5507 = readmatrix('Data to Analyze/5507.cnmfe.LT-imaging-T1.csv');

% get the time of the spikes 
timeOfSpikes5507 = table2array(spikes5507(1:end,1:1));

%vvvv code to get the ez track and inscopix times aligned vvvv
FR=14.255; %FR for specific mouse
inscopixStartTimeInSeconds = 172/FR;
ezTrackStart=470/FR; %in seconds
ezTrackEnd= 16*60+2;
inscopixStartTimesToDelete = ezTrackStart-inscopixStartTimeInSeconds;
inscopixEndTimesToDelete=ezTrackEnd-inscopixStartTimeInSeconds;
endOfVidTimeInSecs = 16*60+2-inscopixStartTimeInSeconds;
indices = find(timeOfSpikes5507<inscopixStartTimesToDelete);
indicesEndOfVid=find(timeOfSpikes5507>inscopixEndTimesToDelete);
indices =[indices;indicesEndOfVid];

spikes5507Corrected = spikes5507;
spikes5507Corrected(indices,:)=[];

%shift all values to be aligned with ez track
spikes5507Corrected(:,1) = spikes5507Corrected(:,1)-inscopixStartTimesToDelete;

spikeTimeCorrected = table2array(spikes5507Corrected(:,1));
lenSpikeTime = height(spikeTimeCorrected);
lenEzTrack = height(ezTrack5507);
xValMatched = zeros(lenSpikeTime,1);

allXVals = table2array(ezTrack5507(:,8));
allTimesEzTrack = table2array(ezTrack5507(:,12));

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

spikesWXVal5507 = [spikes5507Corrected xValMatchedTable];

%writetable(spikesWXVal5507,'spikesWXVal5507.csv');

timeOfSpikesCorrect = table2array(spikesWXVal5507(1:end,1:1));

for INDEX =1:1000
   disp("vvvv Loop # vvvv") 
   disp(INDEX);
% want to add random number to the time values in spikesWXVal5508
randomNum=rand*max(spikesWXVal5507{:,1});
ogMax=max(spikesWXVal5507{:,1});
spikesWXVal5507{:,1}=spikesWXVal5507{:,1}+randomNum;
[rows,cols]=find(spikesWXVal5507{:,1}>ogMax);
spikesWXVal5507{rows,1}=spikesWXVal5507{rows,1}-ogMax;

%plot(allTimesEzTrack,allXVals,'.');
%hold on 
%plot(timeOfSpikesCorrect,xValMatched,'ro');

%Bin the code into 10 horizontal bins
% want to get time each mouse is in bin (add 1 for each frame in specific
% bin
% also want the amount of cell firings for each cell in each bin

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

%now want to attach left or right direction to spike data based on the
%spike point in time
leftOrRight_spikes=strings([1,height(spikesWXVal5507)]); 
%^^ that has L or R for spikeWXVal5507
for i = 1:height(spikesWXVal5507)
    %disp(i);
    %disp("Line 95");
    [row,col]=find(allTimesEzTrack>spikesWXVal5507{i,1}-.5);
    for j = row(1):row(end)
        if (spikesWXVal5507{i,1}-allTimesEzTrack(j))<=0.08
            leftOrRight_spikes(i)=leftOrRight(j);
            break
        end
    end
end
%}
% Code to start the Spatial Map
% this code gets the time the mouse spends in each 160 pixel bin
%should make bin about each 160 secs for 5507 at least
startVal = 0;
endVal = 160;
%timeInBinX = [];
timeInBinXR = [];
timeInBinXL = [];

labelArray="";
%cellNames=unique(spikes5507(:,2)); %array of unique cell names in spike data
%eventsForCellJinBinI = [];
for i = 1:10
    labelArray(i)=strcat(num2str(startVal), ' - ', num2str(endVal));
    timeInBinXR(i)=0;
    timeInBinXL(i)=0;
    for j = 1:length(allXVals)
        %eventsForCellJinBinI(i,j)=0;
        % in here loop thru all the cells and then fill up that array
        if allXVals(j)>=startVal && allXVals(j)<endVal && leftOrRight(j)=='L'
            timeInBinXL(i)=timeInBinXL(i)+1;
        elseif allXVals(j)>=startVal && allXVals(j)<endVal && leftOrRight(j)=='R'
            timeInBinXR(i)=timeInBinXR(i)+1;
        end
    end
    startVal=startVal+160;
    endVal = endVal+ 160;
end
%}
if allXVals(159)>=0 && allXVals(159)<160
    disp("true");
end

%tableWTimeInBins = array2table(timeInBinX,'VariableNames',labelArray); %Has frames in each bin

%% this code will get the number of spikes of cell i in bin j

cellNames=unique(spikes5507(:,2)); %array of unique cell names in spike data
spikesForCellInBin=zeros(height(cellNames),10);

spikesForCellInBinR=zeros(height(cellNames),10);
spikesForCellInBinL=zeros(height(cellNames),10);
for B = 1:height(cellNames)
startVal = 0;
endVal = 160;
%disp(B);
%disp("Line 149");
for i = 1:10
    [rows,cols]=find(startVal<= spikesWXVal5507{:,4} &  spikesWXVal5507{:,4}<endVal & ismember(spikesWXVal5507{:,2},char(cellNames{B,1})) );
    for j=1:length(rows)
        if leftOrRight_spikes(j)=='R'
            spikesForCellInBinR(B,i)=spikesForCellInBinR(B,i)+1;

        else
            spikesForCellInBinL(B,i)=spikesForCellInBinL(B,i)+1;
        end
    end
    startVal = startVal+160;
    endVal = endVal+160;
end
end

%get event rate per bin matrix
mouseSecsInBinXL=timeInBinXL/FR;
mouseSecsInBinXR=timeInBinXR/FR;
eventRatePerBinMatxL=spikesForCellInBinL./mouseSecsInBinXL;
eventRatePerBinMatxR=spikesForCellInBinR./mouseSecsInBinXR;
eventRatePerBinMatx=[eventRatePerBinMatxL,flip(eventRatePerBinMatxR,2)]; %here it is

%code to get probability mouse in ith bin 
for i =1:10
    Pi5507L(i) =mouseSecsInBinXL(i)/(sum(mouseSecsInBinXL));
    Pi5507R(i) =mouseSecsInBinXR(i)/(sum(mouseSecsInBinXR));
end

Pi5507 =[Pi5507L,flip(Pi5507R,2)]; %here it is

infoScore5507=infoScore(eventRatePerBinMatx, Pi5507);
infoScore5507=transpose(infoScore5507);
infoScoreMatx(:,INDEX)=infoScore5507;

end

save("1k Info Scores 5507")