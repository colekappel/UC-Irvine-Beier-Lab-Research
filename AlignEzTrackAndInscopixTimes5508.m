clear

%% Note: timeInBinX array is in frames, have to div by FR to get in seconds

spikes5508 = readtable('Data to Analyze/5508.spikes.LT-imaging-T1.csv');
ezTrack5508 = readtable('Data to Analyze/WIN_20230601_5508_LT_Imaging_T1_LocationOutput.csv');
%cnmfe5508 = readmatrix('Data to Analyze/5508.cnmfe.LT-imaging-T1.csv');

% get the time of the spikes 
timeOfSpikes5508 = table2array(spikes5508(1:end,1:1));

%vvvv code to get the ez track and inscopix times aligned vvvv
FR=14.216;
inscopixStartTimeInSeconds = 55/FR;
ezTrackStart=270/FR;
ezTrackEnd=13804/FR;
inscopixStartTimesToDelete = ezTrackStart-inscopixStartTimeInSeconds;
inscopixEndTimesToDelete=ezTrackEnd-inscopixStartTimeInSeconds;
indices = find(timeOfSpikes5508<inscopixStartTimesToDelete);
endOfVidTimeInSecs = 16*60+11-inscopixStartTimeInSeconds;
indicesEndOfVid=find(timeOfSpikes5508>inscopixEndTimesToDelete);
indices = [indices;indicesEndOfVid];

spikes5508Corrected = spikes5508;
spikes5508Corrected(indices,:)=[];

%shift all values to be aligned with ez track
spikes5508Corrected(:,1) = spikes5508Corrected(:,1)-inscopixStartTimesToDelete;

min(spikes5508Corrected(:,1))

spikeTimeCorrected = table2array(spikes5508Corrected(:,1));
lenSpikeTime = height(spikeTimeCorrected);
lenEzTrack = height(ezTrack5508);
xValMatched = zeros(lenSpikeTime,1);

allXVals = table2array(ezTrack5508(:,8));
allTimesEzTrack = table2array(ezTrack5508(:,12));

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

spikesWXVal5508 = [spikes5508Corrected xValMatchedTable];

%writetable(spikesWXVal5508,'spikesWXVal5508.csv');

%Now want to make a dv/dt array so we can label mouse as moving left or
%right

%this for loop will label xvals starting from 2 to (end - 1) th point
%if its not moving we say its moving left arbitrarily and also at end
%points say its moving left arbitrarily
for INDEX =1:1000
   disp("vvvv Loop # vvvv") 
   disp(INDEX);
% want to add random number to the time values in spikesWXVal5508
randomNum=rand*max(spikesWXVal5508{:,1});
ogMax=max(spikesWXVal5508{:,1});
spikesWXVal5508{:,1}=spikesWXVal5508{:,1}+randomNum;
[rows,cols]=find(spikesWXVal5508{:,1}>ogMax);
spikesWXVal5508{rows,1}=spikesWXVal5508{rows,1}-ogMax;
%{
for i = 1:height(spikesWXVal5508)
    if spikesWXVal5508{i,1}>ogMax
        spikesWXVal5508{i,1}=spikesWXVal5508{i,1}-ogMax;
    end
end
%}
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
leftOrRight_spikes=strings([1,height(spikesWXVal5508)]); 
%^^ that has L or R for spikeWXVal5508
for i = 1:height(spikesWXVal5508)
    %disp(i);
    %disp("Line 95");
    [row,col]=find(allTimesEzTrack>spikesWXVal5508{i,1}-.5);
    for j = row(1):row(end)
        if (spikesWXVal5508{i,1}-allTimesEzTrack(j))<=0.08
            leftOrRight_spikes(i)=leftOrRight(j);
            break
        end
    end
end

%% Code to start the Spatial Map
%% this code gets the time the mouse spends in each 160 pixel bin
%should make bin about each 160 secs for 5508 at least
startVal = 0;
endVal = 160;
%timeInBinX = [];

timeInBinXR = [];
timeInBinXL = [];
labelArray="";
%cellNames=unique(spikes5508(:,2)); %array of unique cell names in spike data
%eventsForCellJinBinI = [];
for i = 1:10
    labelArray(i)=strcat(num2str(startVal), ' - ', num2str(endVal));
    timeInBinXR(i)=0;
    timeInBinXL(i)=0;
    %timeInBinX(i)=0;
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


%tableWTimeInBins = array2table(timeInBinX,'VariableNames',labelArray); %Has frames in each bin

%% this code will get the number of events of cell i in bin j

cellNames=unique(spikes5508(:,2)); %array of unique cell names in spike data
%spikesForCellInBin=zeros(height(cellNames),10);

spikesForCellInBinR=zeros(height(cellNames),10);
spikesForCellInBinL=zeros(height(cellNames),10);
%startVal = 0;
%endVal = 160;
for B = 1:height(cellNames)
startVal = 0;
endVal = 160;
%disp(B);
%disp("Line 141");
for i = 1:10
    [rows,cols]=find(startVal<= spikesWXVal5508{:,4} &  spikesWXVal5508{:,4}<endVal & ismember(spikesWXVal5508{:,2},char(cellNames{B,1})) );
    %[moreRows,moreCols]=find(char(spikesWXVal5508{:,2})==char(cellNames{B,1}));
    for j=1:length(rows)
        if  leftOrRight_spikes(j)=='R'
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
    Pi5508L(i) =mouseSecsInBinXL(i)/(sum(mouseSecsInBinXL));
    Pi5508R(i) =mouseSecsInBinXR(i)/(sum(mouseSecsInBinXR));
end

Pi5508 =[Pi5508L,flip(Pi5508R,2)]; %here it is

infoScore5508=infoScore(eventRatePerBinMatx, Pi5508);
infoScore5508=transpose(infoScore5508);
infoScoreMatx(:,INDEX)=infoScore5508;
end

save("1k Info Scores")