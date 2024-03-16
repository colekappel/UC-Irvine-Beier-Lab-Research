clear 

%5xFAD+
%% cohort 1: 0 (cohort ID)
mouse5507=load('Mouse Data/Data_5507EndBinsRemoved.mat');
%mouse5394=load('Mouse Data/Data_5394EndBinsRemoved.mat');
mouse5511=load('Mouse Data/Data_5511EndBinsRemoved.mat');

%% cohort 2: 1
mouse6025=load('Mouse Data/Data_6025EndBinsRemoved.mat');
mouse6026=load('Mouse Data/Data_6026EndBinsRemoved.mat');
mouse6027=load('Mouse Data/Data_6027EndBinsRemoved.mat');
mouse6028=load('Mouse Data/Data_6028EndBinsRemoved.mat');
mouse6029=load('Mouse Data/Data_6029EndBinsRemoved.mat');

%5xFAD-
%% cohort 1: 0
mouse5508=load('Mouse Data/Data_5508EndBinsRemoved.mat');
mouse5391=load('Mouse Data/Data_5391EndBinsRemoved.mat');
mouse5399=load('Mouse Data/Data_5399EndBinsRemoved.mat');
mouse5455=load('Mouse Data/Data_5455EndBinsRemoved.mat');

%% cohort 2: 1
%mouse6039=load('Mouse Data/Data_6039EndBinsRemoved.mat');
mouse6041=load('Mouse Data/Data_6041EndBinsRemoved.mat');

%lmeFormula = '~FadPlusorMinus+(1|MouseID)+(1|cohortID)';
lmeFormula = '~FadPlusorMinus+(1|MouseID)';

%
%5xFAD- Diagonal Line:
placeCellsArr_5xMinus=[mouse5508.placeCellsArr;mouse5391.placeCellsArr;mouse5399.placeCellsArr;...
    mouse5455.placeCellsArr;mouse6041.placeCellsArr];

%code to visualize diagonal line:
%
diagMinus=figure;
orderMatx_5xMinus=[];

for i =1:height(placeCellsArr_5xMinus)
    [m,j]=max(placeCellsArr_5xMinus(i,:)); %find max in row
    orderMatx_5xMinus(i)=j; %j the column of max val i think
end

[B,I]=sort(orderMatx_5xMinus);
sortedPlaceCells_5xMinus=[];

for i = 1: length(orderMatx_5xMinus)
    sortedPlaceCells_5xMinus(i,:)=placeCellsArr_5xMinus(I(i),:)./max(placeCellsArr_5xMinus(I(i),:));
end
%
%h_minus=heatmap(sortedPlaceCells_5xMinus(35:360,:)); % 5x plus diagonal

h_minus=heatmap(sortedPlaceCells_5xMinus,'Colormap',hot,'ColorLimits',[0 1]);
h_minus.GridVisible='off';
title("5xFAD- Place Cell Visualization")
ylabel("Place Cells")
xlabel("Event Rate Per Bin")
saveas(diagMinus,'MinusDiagonal.svg')
%
diagPlus=figure;
%5xFAD+ Diagonal Line:
placeCellsArr_5xPlus=[mouse5507.placeCellsArr;mouse5511.placeCellsArr;...
    mouse6025.placeCellsArr;mouse6026.placeCellsArr;mouse6027.placeCellsArr;...
    mouse6028.placeCellsArr;mouse6029.placeCellsArr];

%code to visualize diagonal line:
%
orderMatx_5xPlus=[];

for i =1:height(placeCellsArr_5xPlus)
    [m,j]=max(placeCellsArr_5xPlus(i,:)); %find max in row
    orderMatx_5xPlus(i)=j; %j the column of max val i think
end

[B,I]=sort(orderMatx_5xPlus);
sortedPlaceCells_5xPlus=[];
for i = 1: length(orderMatx_5xPlus)
    sortedPlaceCells_5xPlus(i,:)=placeCellsArr_5xPlus(I(i),:)./max(placeCellsArr_5xPlus(I(i),:));
end

%h_plus=heatmap(sortedPlaceCells_5xPlus(62:476,:)); % 5x plus diagonal
h_plus=heatmap(sortedPlaceCells_5xPlus,'Colormap',hot,'ColorLimits',[0 1]);
h_plus.GridVisible='off';
title("5xFAD+ Place Cell Visualization")
ylabel("Place Cells")
xlabel("Event Rate Per Bin")
saveas(diagPlus,'PlusDiagonal.svg')

%}
%{
f5507 = figure;
heatmap(mouse5507.sortedPlaceCells,'Colormap',hot)
ylabel('Cell')
xlabel('Bin (End Bins have been Removed)')
saveas(f5507,'5507Diagonal.svg')

f5508 = figure;
heatmap(mouse5508.sortedPlaceCells,'Colormap',hot)
ylabel('Cell')
xlabel('Bin (End Bins have been Removed)')
saveas(f5508,'5508Diagonal.svg')
figure;
%}
%{
StdAroundTurns= std(mouse5507.top10MeanValsAroundTurns,0,1);
top50MeanAcrossAllTurns = mean(mouse5507.top10MeanValsAroundTurns,1);
timeS = linspace(-1.5,1.5,length(top50MeanAcrossAllTurns));
errorbar(timeS,top50MeanAcrossAllTurns,StdAroundTurns); % plot error bars
%plot(timeS,top50MeanAcrossAllTurns);
title('Average Amplitude of the 50 cells with Highest Event Rate Averaged over all Turns')
xlabel('Time in Seconds (turn happens at t=0)')
ylabel('Amplitude df/f')
%}
%% same code as above but for individual neurons
%{
neuron = mouse5507.seventhNeuronValAtTurns;
StdAroundTurns= std(neuron,0,1);
neuronMeanAcrossAllTurns = mean(neuron,1);
timeS = linspace(-1.5,1.5,length(neuronMeanAcrossAllTurns));
errorbar(timeS,neuronMeanAcrossAllTurns,StdAroundTurns); % plot error bars
%plot(timeS,neuronMeanAcrossAllTurns);
title('Average Amplitude of a Single Neuron with a High Event Rate Averaged over all Turns')
xlabel('Time in Seconds (turn happens at t=0)')
ylabel('Amplitude df/f')
%}
%{
%% to plot 5507 amplitude vs turns
plot(mouse5507.timesAroundTurns(14,:),mouse5507.top10MeanValsAroundTurns(14,:));
title('Average Amplitude of the 50 cells with Highest Event Rate at Turns: 14th Turn')
xlabel('Time in Seconds (turn happens in the middle of the plot)')
ylabel('Amplitude df/f')
figure;

plot(mouse5507.timesAroundTurns(14,:),mouse5507.top10MeanValsAroundTurns(10,:));
title('Average Amplitude of the 50 cells with Highest Event Rate at Turns: 10th Turn')
xlabel('Time in Seconds (turn happens in the middle of the plot)')
ylabel('Amplitude df/f')
figure;

plot(mouse5507.timesAroundTurns(14,:),mouse5507.top10MeanValsAroundTurns(12,:));
title('Average Amplitude of the 50 cells with Highest Event Rate at Turns: 12th Turn')
xlabel('Time in Seconds (turn happens in the middle of the plot)')
ylabel('Amplitude df/f')
%}

%% for time per crossing stats

%}
%{
%% bar charts for stability
%% bar chart for stability Average
bar([mean(allStabilityAvMinus) 0])
hold on
bar([0 mean(allStabilityAvPlus)])
xlabel("Condition",'FontSize',20);
ylabel("Stability",'FontSize',20);
title("Bar chart for Stability: Average",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

%% bar chart for stability Odd vs Even
figure;
bar([mean(allStabilityOEMinus) 0])
hold on
bar([0 mean(allStabilityOEPlus)])
xlabel("Condition",'FontSize',20);
ylabel("Stability",'FontSize',20);
title("Bar chart for Stability: Odd vs Even Trials",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

%% bar chart for stability 1v2
figure;
bar([mean(allStability1v2Minus) 0])
hold on
bar([0 mean(allStability1v2Plus)])
xlabel("Condition",'FontSize',20);
ylabel("Stability",'FontSize',20);
title("Bar chart for Stability: 1st vs 2nd half",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')
%}

%}
%
%% p tests, ks tests and LME for info score, sparsity, spatial coherence, stability
%% amplitude and event rate
%
%% For stability average
allStabilityAvPlus = [mouse5507.stabilityAverage;mouse5511.stabilityAverage; mouse6025.stabilityAverage;...
    mouse6026.stabilityAverage ; mouse6027.stabilityAverage ; mouse6028.stabilityAverage; mouse6029.stabilityAverage ];
allStabilityAvMinus = [mouse5508.stabilityAverage;mouse5391.stabilityAverage;mouse5399.stabilityAverage;mouse5455.stabilityAverage;...
     mouse6041.stabilityAverage];
[hStabAvKs,pStabAvKs] = kstest2(allStabilityAvPlus,allStabilityAvMinus);
[hStabAvTt,pStabAvTt] = ttest2(allStabilityAvPlus,allStabilityAvMinus);
meanStabAvplus = mean(allStabilityAvPlus);
meanStabAvminus = mean(allStabilityAvMinus);

%% For stability 1v2
allStability1v2Plus = [mouse5507.stability1v2;mouse5511.stability1v2; mouse6025.stability1v2;...
    mouse6026.stability1v2 ; mouse6027.stability1v2 ; mouse6028.stability1v2; mouse6029.stability1v2];
allStability1v2Minus = [mouse5508.stability1v2;mouse5391.stability1v2;mouse5399.stability1v2;mouse5455.stability1v2;...
     mouse6041.stability1v2];
[hStab1v2Ks,pStab1v2Ks] = kstest2(allStability1v2Plus,allStability1v2Minus);
[hStab1v2Tt,pStab1v2Tt] = ttest2(allStability1v2Plus,allStability1v2Minus);
meanStab1v2plus = mean(allStability1v2Plus);
meanStab1v2minus = mean(allStability1v2Minus);

%% For stability odd v even
allStabilityOEPlus = [mouse5507.stabilityOddEven;mouse5511.stabilityOddEven; mouse6025.stabilityOddEven;...
    mouse6026.stabilityOddEven ; mouse6027.stabilityOddEven ; mouse6028.stabilityOddEven; mouse6029.stabilityOddEven];
allStabilityOEMinus = [mouse5508.stabilityOddEven;mouse5391.stabilityOddEven;mouse5399.stabilityOddEven;mouse5455.stabilityOddEven;...
     mouse6041.stabilityOddEven];
[hStabOEKs,pStabOEKs] = kstest2(allStabilityOEPlus,allStabilityOEMinus);
[hStabOETt,pStabOETt] = ttest2(allStabilityOEPlus,allStabilityOEMinus);
meanStabOEplus = mean(allStabilityOEPlus);
meanStabOEminus = mean(allStabilityOEMinus);

%% info score: .theInfoScore
%
allInfoScoresPlus = [mouse5507.theInfoScore ; mouse5511.theInfoScore; mouse6025.theInfoScore;...
    mouse6026.theInfoScore ; mouse6027.theInfoScore ; mouse6028.theInfoScore ; mouse6029.theInfoScore];
allInfoScoresMinus = [mouse5508.theInfoScore; mouse5391.theInfoScore; mouse5399.theInfoScore; mouse5455.theInfoScore;...
     mouse6041.theInfoScore];
[hInfoKs,pInfoKs] = kstest2(allInfoScoresPlus,allInfoScoresMinus);
[hInfoTt,pInfoTt] = ttest2(allInfoScoresPlus,allInfoScoresMinus); %%% here on editing
%}

%% sparsity: .sparsity
allSparsPlus = [mouse5507.sparsity;mouse5511.sparsity; mouse6025.sparsity;...
    mouse6026.sparsity ; mouse6027.sparsity ; mouse6028.sparsity; mouse6029.sparsity];
allSparssMinus  = [mouse5508.sparsity;mouse5391.sparsity;mouse5399.sparsity;mouse5455.sparsity;...
     mouse6041.sparsity];
[hSparsKs,pSparsKs] = kstest2(allSparsPlus,allSparssMinus);
[hSparsTt,pSparsTt] = ttest2(allSparsPlus,allSparssMinus);

%% spatial coherence: .spatialCoherence
%muting bc we dont use it
%{ 
allSpatCPlus = [mouse5507.spatialCoherence;mouse5394.spatialCoherence;mouse5511.spatialCoherence; mouse6025.spatialCoherence;...
    mouse6026.spatialCoherence ; mouse6027.spatialCoherence ; mouse6028.spatialCoherence; mouse6029.spatialCoherence];
allSpatCMinus  = [mouse5508.spatialCoherence;mouse5391.spatialCoherence;mouse5399.spatialCoherence;mouse5455.spatialCoherence;...
     mouse6041.spatialCoherence];
[hSpatCKs,pSpatCKs] = kstest2(allSpatCPlus,allSpatCMinus);
[hSpatCTt,pSpatCTt] = ttest2(allSpatCPlus,allSpatCMinus);
%}

%% amplitude: .avAmplitudeForEachCell
 PlusAmpAllCells = [mouse5507.avAmplitudeForEachCell;mouse5511.avAmplitudeForEachCell; mouse6025.avAmplitudeForEachCell;...
    mouse6026.avAmplitudeForEachCell ; mouse6027.avAmplitudeForEachCell ; mouse6028.avAmplitudeForEachCell; mouse6029.avAmplitudeForEachCell]; 
 MinusAmpAllCells  = [mouse5508.avAmplitudeForEachCell;mouse5391.avAmplitudeForEachCell;mouse5399.avAmplitudeForEachCell;mouse5455.avAmplitudeForEachCell;...
     mouse6041.avAmplitudeForEachCell]; 
 [hAmpKs, pAmpKs ] = kstest2(MinusAmpAllCells,PlusAmpAllCells) ;
 [hAmpTt, pAmpTt ] = ttest2(MinusAmpAllCells,PlusAmpAllCells) ;

 %% event rate: .eventRatePerCell
 PlusERAllCells = [mouse5507.eventRatePerCell;mouse5511.eventRatePerCell; mouse6025.eventRatePerCell;...
    mouse6026.eventRatePerCell ; mouse6027.eventRatePerCell ; mouse6028.eventRatePerCell; mouse6029.eventRatePerCell]; 
 MinusERAllCells  = [mouse5508.eventRatePerCell;mouse5391.eventRatePerCell;mouse5399.eventRatePerCell;mouse5455.eventRatePerCell;...
     mouse6041.eventRatePerCell];  
 [hERKs, pERKs ] = kstest2(MinusERAllCells,PlusERAllCells) ;
 [hERTt, pERTt ] = ttest2(MinusERAllCells,PlusERAllCells) ;
%

%
%% time per crossing LME Model
%concatenate the LME tables into one table 
LMETableTimePerCrossing = [mouse5507.LMETableTimePerCrossing; mouse5511.LMETableTimePerCrossing;...
    mouse6025.LMETableTimePerCrossing; mouse6026.LMETableTimePerCrossing; mouse6027.LMETableTimePerCrossing; mouse6028.LMETableTimePerCrossing;... 
    mouse6029.LMETableTimePerCrossing;  mouse6041.LMETableTimePerCrossing;... 
    mouse5508.LMETableTimePerCrossing; mouse5391.LMETableTimePerCrossing;mouse5399.LMETableTimePerCrossing;mouse5455.LMETableTimePerCrossing];
%LMETable.FadPlusorMinus=string(LMETable{:,2});
%
formulaTPC = strcat('timePerCrossing',lmeFormula);

lmeTPC = fitlme(LMETableTimePerCrossing,formulaTPC);

allTimePerCrossingPlus = [mouse5507.timePerCrossing;mouse5511.timePerCrossing; mouse6025.timePerCrossing;...
    mouse6026.timePerCrossing;mouse6027.timePerCrossing;mouse6028.timePerCrossing; mouse6029.timePerCrossing];
allTimePerCrossingMinus = [mouse5508.timePerCrossing;mouse5391.timePerCrossing;mouse5399.timePerCrossing;mouse5455.timePerCrossing;...
    mouse6041.timePerCrossing];
[hTPCKs,pTPCKs] = kstest2(allTimePerCrossingPlus,allTimePerCrossingMinus);
[hTPCTt,pTPCTt] = ttest2(allTimePerCrossingPlus,allTimePerCrossingMinus);

%% Stability ks test, t test and lme test p values:
%

%% LME Models for stability:
%% LME Model for stability 1v2
%
LMETableStability1v2 = [mouse5507.LMETableStability1v2; mouse5511.LMETableStability1v2;...
    mouse6025.LMETableStability1v2; mouse6026.LMETableStability1v2; mouse6027.LMETableStability1v2; mouse6028.LMETableStability1v2;...
    mouse6029.LMETableStability1v2;  mouse6041.LMETableStability1v2;...
    mouse5508.LMETableStability1v2; mouse5391.LMETableStability1v2;mouse5399.LMETableStability1v2;mouse5455.LMETableStability1v2];
%LMETable.FadPlusorMinus=string(LMETable{:,2});

formulaStability1v2 = strcat('stability1v2',lmeFormula);
%% new lme: 'stability1v2~FadPlusorMinus+(1|MouseID)+(1|CohortID)'
lmeStability1v2= fitlme(LMETableStability1v2,formulaStability1v2);

%% LME Model for stability OE: LMETableStabilityOE
LMETableStabilityOE = [mouse5507.LMETableStabilityOE; mouse5511.LMETableStabilityOE;...
    mouse6025.LMETableStabilityOE; mouse6026.LMETableStabilityOE; mouse6027.LMETableStabilityOE; mouse6028.LMETableStabilityOE;...
    mouse6029.LMETableStabilityOE;  mouse6041.LMETableStabilityOE;...
    mouse5508.LMETableStabilityOE; mouse5391.LMETableStabilityOE;mouse5399.LMETableStabilityOE;mouse5455.LMETableStabilityOE];
%LMETable.FadPlusorMinus=string(LMETable{:,2});

formulaStabilityOE = strcat('stabilityOddEven',lmeFormula);

lmeStabilityOE= fitlme(LMETableStabilityOE,formulaStabilityOE);

%% LME Model for stability Average: LMETableStabilityAv
LMETableStabilityAv = [mouse5507.LMETableStabilityAv; mouse5511.LMETableStabilityAv;...
    mouse6025.LMETableStabilityAv; mouse6026.LMETableStabilityAv; mouse6027.LMETableStabilityAv; mouse6028.LMETableStabilityAv;...
    mouse6029.LMETableStabilityAv;  mouse6041.LMETableStabilityAv;...
    mouse5508.LMETableStabilityAv; mouse5391.LMETableStabilityAv;mouse5399.LMETableStabilityAv;mouse5455.LMETableStabilityAv];
%LMETable.FadPlusorMinus=string(LMETable{:,2});

formulaStabilityAv = strcat('stabilityAverage',lmeFormula);

lmeStabilityAv= fitlme(LMETableStabilityAv,formulaStabilityAv);

%% Amplitude LME Model: .LMETableAmplitude
%concatenate the LME tables into one table 
LMETableAmplitude= [mouse5507.LMETableAmplitude; mouse5511.LMETableAmplitude;...
    mouse6025.LMETableAmplitude; mouse6026.LMETableAmplitude; mouse6027.LMETableAmplitude; mouse6028.LMETableAmplitude;...
    mouse6029.LMETableAmplitude; mouse6041.LMETableAmplitude;...
    mouse5508.LMETableAmplitude; mouse5391.LMETableAmplitude;mouse5399.LMETableAmplitude;mouse5455.LMETableAmplitude];

%LMETableAmplitude.FadPlusorMinus=categorical(LMETableAmplitude{:,2});
%
formulaAmp = strcat('avAmplitudeForEachCell',lmeFormula);

lmeAmp = fitlme(LMETableAmplitude,formulaAmp);

%% Event Rate LME Model: .LMETableEventRate

LMETableER = [mouse5507.LMETableEventRate; mouse5511.LMETableEventRate;...
    mouse6025.LMETableEventRate; mouse6026.LMETableEventRate; mouse6027.LMETableEventRate; mouse6028.LMETableEventRate;...
    mouse6029.LMETableEventRate; mouse6041.LMETableEventRate;...
    mouse5508.LMETableEventRate; mouse5391.LMETableEventRate;mouse5399.LMETableEventRate;mouse5455.LMETableEventRate];
%LMETable.FadPlusorMinus=string(LMETable{:,2});

formulaER = strcat('eventRatePerCell',lmeFormula);

lmeER= fitlme(LMETableER,formulaER);

%% Info Score LME Model: .LMETableInfoScore

LMETableInfo  = [mouse5507.LMETableInfoScore; mouse5511.LMETableInfoScore;...
    mouse6025.LMETableInfoScore; mouse6026.LMETableInfoScore; mouse6027.LMETableInfoScore; mouse6028.LMETableInfoScore;...
    mouse6029.LMETableInfoScore;  mouse6041.LMETableInfoScore;...
    mouse5508.LMETableInfoScore; mouse5391.LMETableInfoScore;mouse5399.LMETableInfoScore;mouse5455.LMETableInfoScore];
%LMETable.FadPlusorMinus=string(LMETable{:,2});

formulaInfo = strcat('theInfoScore',lmeFormula);

lmeInfo= fitlme(LMETableInfo,formulaInfo);

%% Sparsity LME Model: .LMETableSparsity

LMETableSparsity  = [mouse5507.LMETableSparsity; mouse5511.LMETableSparsity;...
    mouse6025.LMETableSparsity; mouse6026.LMETableSparsity; mouse6027.LMETableSparsity; mouse6028.LMETableSparsity;...
    mouse6029.LMETableSparsity; mouse6041.LMETableSparsity;...
    mouse5508.LMETableSparsity; mouse5391.LMETableSparsity;mouse5399.LMETableSparsity;mouse5455.LMETableSparsity];
%LMETable.FadPlusorMinus=string(LMETable{:,2});

formulaSparsity = strcat('sparsity',lmeFormula);

lmeSparsity= fitlme(LMETableSparsity,formulaSparsity);

%% Sparial Coherence LME Model: .LMETableSpatC
%mute since we dont use it
%{
LMETableSpatC = [mouse5507.LMETableSpatC; mouse5394.LMETableSpatC;mouse5511.LMETableSpatC;...
    mouse6025.LMETableSpatC; mouse6026.LMETableSpatC; mouse6027.LMETableSpatC; mouse6028.LMETableSpatC;...
    mouse6029.LMETableSpatC;  mouse6041.LMETableSpatC;...
    mouse5508.LMETableSpatC; mouse5391.LMETableSpatC;mouse5399.LMETableSpatC;mouse5455.LMETableSpatC];
%LMETable.FadPlusorMinus=string(LMETable{:,2});

formulaSpatC = strcat('spatialCoherence',lmeFormula);

lmeSpatC= fitlme(LMETableSpatC,formulaSpatC);
%}
%% should also do ks test and t test for amp and event rate across all cells
%

 %}
%% bar charts for all analyses
%{
bar([mean(allInfoScoresMinus) 0])
hold on
bar([0 mean(allInfoScoresPlus)])
xlabel("Condition",'FontSize',20);
ylabel("Information Score",'FontSize',20);
title("Bar chart for Average Information Score",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;
bar([mean(MinusERAllCells) 0])
hold on
bar([0 mean(PlusERAllCells)])
xlabel("Condition",'FontSize',20);
ylabel("Event Rate",'FontSize',20);
title("Bar chart for Average Event Rate",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;
bar([mean(allSparssMinus) 0])
hold on
bar([0 mean(allSparsPlus)])
xlabel("Condition",'FontSize',20);
ylabel("Event Rate",'FontSize',20);
title("Bar chart for Average Sparsity",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;
bar([mean(allSpatCMinus) 0])
hold on
bar([0 mean(allSpatCPlus)])
xlabel("Condition",'FontSize',20);
ylabel("Event Rate",'FontSize',20);
title("Bar chart for Average Spatial Coherence",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;
bar([mean(MinusAmpAllCells) 0])
hold on
bar([0 mean(PlusAmpAllCells)])
xlabel("Condition",'FontSize',20);
ylabel("Event Rate",'FontSize',20);
title("Bar chart for Average Amplitude",'FontSize',20)
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')
%}
%
%% cdf plots for sparsity, spatial coherence, info score, amplitude, event
%% rate
%{
figure;
cdfplot(allSparssMinus);
hold on
cdfplot(allSparsPlus);
title("Cumulative Distribution Plot for Sparsity",'FontSize',14)
xlabel("Sparsity",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

%{
figure;
cdfplot(allSpatCMinus);
hold on
cdfplot(allSpatCPlus);
title("Cumulative Distribution Plot for Spatial Coherence",'FontSize',14)
xlabel("Spatial Coherence",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')
%}

figure;
cdfplot(allInfoScoresMinus);
hold on
cdfplot(allInfoScoresPlus);
title("Cumulative Distribution Plot for Information Score",'FontSize',14)
xlabel("Information Score",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;

cdfplot(MinusAmpAllCells);
hold on
cdfplot(PlusAmpAllCells);
title("Cumulative Distribution Plot for Amplitude",'FontSize',14)
xlabel("Average Amplitude Per Cell",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;

cdfplot(MinusERAllCells);
hold on
cdfplot(PlusERAllCells);
title("Cumulative Distribution Plot for Event Rate",'FontSize',14)
xlabel("Event Rate Per Cell",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;
%
%
% cdf plots for stability
cdfplot(allStabilityAvMinus);
hold on
cdfplot(allStabilityAvPlus);
title("Cumulative Distribution Plot for Stability: Average",'FontSize',14)
xlabel("Stability",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;
cdfplot(allStability1v2Minus);
hold on
cdfplot(allStability1v2Plus);
title("Cumulative Distribution Plot for Stability: 1st vs 2nd half",'FontSize',14)
xlabel("Stability",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')

figure;
cdfplot(allStabilityOEMinus);
hold on
cdfplot(allStabilityOEPlus);
title("Cumulative Distribution Plot for Stability: Odd vs Even trials",'FontSize',14)
xlabel("Stability",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')
%}
%{
figure;
cdfplot(allTimePerCrossingMinus);
hold on
cdfplot(allTimePerCrossingPlus);
title("Cumulative Distribution Plot for Time Per Crossing",'FontSize',14)
xlabel("Time Per Crossing",'FontSize',14);
ylabel("% of cells",'FontSize',14);
legend('5xFAD-','5xFAD+','FontSize',20,'Location','best')
%}

%% box plots for event rate, amplitude, info score, sparsity  

%% box plots for event rate
%
figure;
gg1 = repmat({'5xFAD-'},length(MinusERAllCells),1);
gg2= repmat({'5xFAD+'},length(PlusERAllCells),1);
gg=[gg1;gg2];
xx=[MinusERAllCells;PlusERAllCells];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Event Rate Per Cell for all Cells Across all Mice Comparison','FontSize',20)
ylabel('Event Rate','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
figure;
%% box plot for amplitude
gg1 = repmat({'5xFAD-'},length(MinusAmpAllCells),1);
gg2= repmat({'5xFAD+'},length(PlusAmpAllCells),1);
gg=[gg1;gg2];
xx=[MinusAmpAllCells;PlusAmpAllCells];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Average Amplitude Per Cell for all Cells Across all Mice Comparison','FontSize',20)
ylabel('Average Amplitude','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
figure;
%% box plot for info score
gg1 = repmat({'5xFAD-'},length(allInfoScoresMinus),1);
gg2= repmat({'5xFAD+'},length(allInfoScoresPlus),1);
gg=[gg1;gg2];
xx=[allInfoScoresMinus ; allInfoScoresPlus];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Information Score for all Cells Across all Mice Comparison','FontSize',20)
ylabel('Information Score','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
figure;
%% box plot for sparsity
gg1 = repmat({'5xFAD-'},length(allSparssMinus),1);
gg2= repmat({'5xFAD+'},length(allSparsPlus),1);
gg=[gg1;gg2];
xx=[allSparssMinus;allSparsPlus];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Sparsity for all Cells Across all Mice Comparison','FontSize',20)
ylabel('Sparsity','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
figure;

%% box plot for spatial coherence
%{
gg1 = repmat({'5xFAD-'},length(allSpatCMinus),1);
gg2= repmat({'5xFAD+'},length(allSpatCPlus),1);
gg=[gg1;gg2];
xx=[allSpatCMinus;allSpatCPlus];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Spatial Coherence for all Cells Across all Mice Comparison','FontSize',20)
ylabel('Spatial Coherence','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

figure
%}
%% box plots for stabilty
%
%% box plot for stabilty Average
gg1 = repmat({'5xFAD-'},length(allStabilityAvMinus),1);
gg2= repmat({'5xFAD+'},length(allStabilityAvPlus),1);
gg=[gg1;gg2];
xx=[allStabilityAvMinus;allStabilityAvPlus];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Box Plot for Stability: Average','FontSize',20)
ylabel('Stability','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
figure;
%% box plot for stabilty 1 v 2
gg1 = repmat({'5xFAD-'},length(allStability1v2Minus),1);
gg2= repmat({'5xFAD+'},length(allStability1v2Plus),1);
gg=[gg1;gg2];
xx=[allStability1v2Minus;allStability1v2Plus];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Box Plot for Stability: 1st vs 2nd half','FontSize',20)
ylabel('Stability','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
figure;
%% box plot for stabilty Odd v Even
gg1 = repmat({'5xFAD-'},length(allStabilityOEMinus),1);
gg2= repmat({'5xFAD+'},length(allStabilityOEPlus),1);
gg=[gg1;gg2];
xx=[allStabilityOEMinus;allStabilityOEPlus];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Box Plot for Stability: Odd vs Even Trials','FontSize',20)
ylabel('Stability','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
%}
%{
%% box plot for Time per Crossing
gg1 = repmat({'5xFAD-'},length(allTimePerCrossingMinus),1);
gg2= repmat({'5xFAD+'},length(allTimePerCrossingPlus),1);
gg=[gg1;gg2];
xx=[allTimePerCrossingMinus;allTimePerCrossingPlus];
boxplot(xx,gg);
%[hSparsity,pSparsity] = ttest2(mouse5507.sparsity, mouse5508.sparsity);
title('Box Plot for Time Per Crossing','FontSize',20)
ylabel('Time Per Crossing','FontSize',14)

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [red;blue];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

%}
%}
%save("LMEDataEndBinsRemoved")
%}


%writetable(LMETableTimePerCrossing,'/Users/colekappel/Desktop/Best Day Time Per Crossing Stats/TimePerCrossing.csv')

%{
writetable(LMETableAmplitude,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/AmplitudeData.csv')
writetable(LMETableInfo,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/InfoScoreData.csv')
writetable(LMETableSparsity,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/SparsityData.csv')
writetable(LMETableER,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/EventRateData.csv')

writetable(LMETableStability1v2,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/Stability1v2Data.csv')
writetable(LMETableStabilityOE,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/StabilityOE.csv')
writetable(LMETableStabilityAv,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/StabilityAvData.csv')

%writematrix(sortedPlaceCells_5xMinus,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/Diagonal_5xMinus.csv')
%writematrix(sortedPlaceCells_5xPlus,'/Users/colekappel/Desktop/Best Day PC w Excl - Correct/Diagonal_5xPlus.csv')
%^^ heat map csv's

%writetable(LMETableSpatC,'/Users/colekappel/Desktop/Figures and Data for Balbina/SpatCohData.csv')
%}

%{
writematrix(transpose(mouse5391.timePerCrossing),'/Users/colekappel/Desktop/Behavior/5391TimePerCrossing.csv')
writematrix(transpose(mouse5394.timePerCrossing),'/Users/colekappel/Desktop/Behavior/5394TimePerCrossing.csv')
writematrix(transpose(mouse5399.timePerCrossing),'/Users/colekappel/Desktop/Behavior/5399TimePerCrossing.csv')
writematrix(transpose(mouse5455.timePerCrossing),'/Users/colekappel/Desktop/Behavior/5455TimePerCrossing.csv')
writematrix(mouse5507.timePerCrossing,'/Users/colekappel/Desktop/Behavior/5507TimePerCrossing.csv')
writematrix(transpose(mouse5508.timePerCrossing),'/Users/colekappel/Desktop/Behavior/5508TimePerCrossing.csv')
writematrix(transpose(mouse5511.timePerCrossing),'/Users/colekappel/Desktop/Behavior/5511TimePerCrossing.csv')
%}