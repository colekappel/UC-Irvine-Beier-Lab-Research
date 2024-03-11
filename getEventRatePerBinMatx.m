function [Pi,mouseSecsInBinX,eventPerBinMatx,eventRatePerBinMatx]...
    = getEventRatePerBinMatx(NUMBEROFBINS,BINSIZE,allXValsWDirection,cellNames,spikesWXVal,leftOrRight_spikes,FR)

%% Code to start the Spatial Map
%% this code gets the time the mouse spends in each 160 pixel bin
startVal = 0;
%BINSIZE = max(spikesWXVal{:,4})/10;
endVal=BINSIZE;
%timeInBinX = [];

timeInBinXR = [];
timeInBinXL = [];
labelArray="";
%cellNames=unique(spikes5508(:,2)); %array of unique cell names in spike data
%eventsForCellJinBinI = [];
for i = 1:NUMBEROFBINS
    labelArray(i)=strcat(num2str(startVal), ' - ', num2str(endVal));
    timeInBinXR(i)=0;
    timeInBinXL(i)=0;
    %timeInBinX(i)=0;
    for j = 1:length(allXValsWDirection)
        %eventsForCellJinBinI(i,j)=0;
        % in here loop thru all the cells and then fill up that array
       if double(allXValsWDirection(j,1))>=startVal && double(allXValsWDirection(j,1))<endVal && allXValsWDirection(j,2)=='L'
            timeInBinXL(i)=timeInBinXL(i)+1;
        elseif double(allXValsWDirection(j,1))>=startVal && double(allXValsWDirection(j,1))<endVal && allXValsWDirection(j,2)=='R'
            timeInBinXR(i)=timeInBinXR(i)+1;
        end
    end
    startVal=startVal+BINSIZE;
    endVal = endVal+ BINSIZE;
end
%}

%tableWTimeInBins = array2table(timeInBinX,'VariableNames',labelArray); %Has frames in each bin

%% this code will get the number of events of cell i in bin j
%spikesForCellInBin=zeros(height(cellNames),10);

spikesForCellInBinR=zeros(height(cellNames),NUMBEROFBINS);
spikesForCellInBinL=zeros(height(cellNames),NUMBEROFBINS);
%startVal = 0;
%endVal = 160;
for B = 1:height(cellNames)
startVal = 0;
endVal = BINSIZE;
%disp(B);
%disp("Line 141");
for i = 1:NUMBEROFBINS
    [rows,cols]=find(startVal<= spikesWXVal{:,4} &  spikesWXVal{:,4}<endVal & ismember(spikesWXVal{:,2},char(cellNames{B,1})) );
    %[moreRows,moreCols]=find(char(spikesWXVal5508{:,2})==char(cellNames{B,1}));
    for j=1:length(rows)
        if  leftOrRight_spikes(rows(j))=='R'
            spikesForCellInBinR(B,i)=spikesForCellInBinR(B,i)+1;

        else
            spikesForCellInBinL(B,i)=spikesForCellInBinL(B,i)+1;
        end
    end
    startVal = startVal+BINSIZE;
    endVal = endVal+BINSIZE;
end
end
%get event rate per bin matrix
mouseSecsInBinXL=timeInBinXL/FR; 
mouseSecsInBinXR=timeInBinXR/FR; 
mouseSecsInBinX=[mouseSecsInBinXL,flip(mouseSecsInBinXR,2)]; %output from func
eventRatePerBinMatxL=spikesForCellInBinL./mouseSecsInBinXL;
eventRatePerBinMatxR=spikesForCellInBinR./mouseSecsInBinXR;
eventRatePerBinMatx=[eventRatePerBinMatxL,flip(eventRatePerBinMatxR,2)]; %output from func
eventPerBinMatx=[spikesForCellInBinL,flip(spikesForCellInBinR,2)];

%code to get probability mouse in ith bin 
for i =1:NUMBEROFBINS
    PiL(i) =mouseSecsInBinXL(i)/(sum(mouseSecsInBinXL));
    PiR(i) =mouseSecsInBinXR(i)/(sum(mouseSecsInBinXR));
end

Pi =[PiL,flip(PiR,2)]; %output from func

Pi(:,[1:NUMBEROFBINS/10, (NUMBEROFBINS+1-NUMBEROFBINS/10):(NUMBEROFBINS+NUMBEROFBINS/10), (NUMBEROFBINS*2+1-NUMBEROFBINS/10):(NUMBEROFBINS*2)])=[];
eventRatePerBinMatx(:,[1:NUMBEROFBINS/10, (NUMBEROFBINS+1-NUMBEROFBINS/10):(NUMBEROFBINS+NUMBEROFBINS/10), (NUMBEROFBINS*2+1-NUMBEROFBINS/10):(NUMBEROFBINS*2)])=[];
eventPerBinMatx(:,[1:NUMBEROFBINS/10, (NUMBEROFBINS+1-NUMBEROFBINS/10):(NUMBEROFBINS+NUMBEROFBINS/10), (NUMBEROFBINS*2+1-NUMBEROFBINS/10):(NUMBEROFBINS*2)])=[];
mouseSecsInBinX(:,[1:NUMBEROFBINS/10, (NUMBEROFBINS+1-NUMBEROFBINS/10):(NUMBEROFBINS+NUMBEROFBINS/10), (NUMBEROFBINS*2+1-NUMBEROFBINS/10):(NUMBEROFBINS*2)])=[];

end
