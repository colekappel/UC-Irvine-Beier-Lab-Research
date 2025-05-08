clear 
%% 5xFAD+
%% cohort 1: 0 (cohort ID)
mouse5507=load('Mouse Data/Data_5507EndBinsRemoved.mat').normedPCArr;
%mouse5394=load('Mouse Data/Data_5394EndBinsRemoved.mat').normedPCArr; left
%out
mouse5511=load('Mouse Data/Data_5511EndBinsRemoved.mat').normedPCArr;
%% cohort 2: 1
mouse6025=load('Mouse Data/Data_6025EndBinsRemoved.mat').normedPCArr;
mouse6026=load('Mouse Data/Data_6026EndBinsRemoved.mat').normedPCArr;
mouse6027=load('Mouse Data/Data_6027EndBinsRemoved.mat').normedPCArr;
mouse6028=load('Mouse Data/Data_6028EndBinsRemoved.mat').normedPCArr;
mouse6029=load('Mouse Data/Data_6029EndBinsRemoved.mat').normedPCArr;
%% 5xFAD-
%% cohort 1: 0
mouse5508=load('Mouse Data/Data_5508EndBinsRemoved.mat').normedPCArr;
mouse5391=load('Mouse Data/Data_5391EndBinsRemoved.mat').normedPCArr;
mouse5399=load('Mouse Data/Data_5399EndBinsRemoved.mat').normedPCArr;
mouse5455=load('Mouse Data/Data_5455EndBinsRemoved.mat').normedPCArr;
%% cohort 2: 1
mouse6039=load('Mouse Data/Data_6039EndBinsRemoved.mat').normedPCArr;
mouse6041=load('Mouse Data/Data_6041EndBinsRemoved.mat').normedPCArr;

pcArrAllMice = [mouse5507;mouse5511;mouse6025;mouse6026;mouse6027;...
    mouse6028;mouse6029;mouse5508;mouse5391;mouse5399;mouse5455;...
    mouse6039;mouse6041];

boolAllMice= pcArrAllMice>0.3;

numPlaceFieldsArrAllMice = zeros(height(boolAllMice),1);
maxWidthArrAllMice = zeros(height(boolAllMice),1);

for J = 1:height(boolAllMice)

bwAllMice=bwareafilt(boolAllMice(J,:),[1,16]);

ccAllMice=bwconncomp(bwAllMice);

numPlaceFieldsAllMice = length(ccAllMice.PixelIdxList); % need this for each cell

numPlaceFieldsArrAllMice(J) = numPlaceFieldsAllMice;

widths5391 = zeros(1,numPlaceFieldsAllMice);
for i = 1:numPlaceFieldsAllMice

widths5391(i)= length(ccAllMice.PixelIdxList{i});

end

maxWidth5391 = max(widths5391); % need this for each cell

maxWidthArrAllMice(J) = maxWidth5391;
end

%% 5xFAD+
%% cohort 1: 0 (cohort ID)
mouse5507lme=load('Mouse Data/Data_5507EndBinsRemoved.mat').LMETableEventRate;
%mouse5394=load('Mouse Data/Data_5394EndBinsRemoved.mat').LMETableEventRate; left
%out
mouse5511lme=load('Mouse Data/Data_5511EndBinsRemoved.mat').LMETableEventRate;
%% cohort 2: 1
mouse6025lme=load('Mouse Data/Data_6025EndBinsRemoved.mat').LMETableEventRate;
mouse6026lme=load('Mouse Data/Data_6026EndBinsRemoved.mat').LMETableEventRate;
mouse6027lme=load('Mouse Data/Data_6027EndBinsRemoved.mat').LMETableEventRate;
mouse6028lme=load('Mouse Data/Data_6028EndBinsRemoved.mat').LMETableEventRate;
mouse6029lme=load('Mouse Data/Data_6029EndBinsRemoved.mat').LMETableEventRate;
%% 5xFAD-
%% cohort 1: 0
mouse5508lme=load('Mouse Data/Data_5508EndBinsRemoved.mat').LMETableEventRate;
mouse5391lme=load('Mouse Data/Data_5391EndBinsRemoved.mat').LMETableEventRate;
mouse5399lme=load('Mouse Data/Data_5399EndBinsRemoved.mat').LMETableEventRate;
mouse5455lme=load('Mouse Data/Data_5455EndBinsRemoved.mat').LMETableEventRate;
%% cohort 2: 1
mouse6039lme=load('Mouse Data/Data_6039EndBinsRemoved.mat').LMETableEventRate;
mouse6041lme=load('Mouse Data/Data_6041EndBinsRemoved.mat').LMETableEventRate;

lmeAllMice = [mouse5507lme;mouse5511lme;mouse6025lme;mouse6026lme;mouse6027lme;...
    mouse6028lme;mouse6029lme;mouse5508lme;mouse5391lme;mouse5399lme;mouse5455lme;...
    mouse6039lme;mouse6041lme];

lmeAllMice.numPlaceFields = numPlaceFieldsArrAllMice;
lmeAllMice.maxPlaceFieldWidth = maxWidthArrAllMice;

lmeAllMice(:,1) = [];

formulaNumPlaceFields =...
    'numPlaceFields~FadPlusorMinus+(1|MouseID)';

lmeNumPlaceFields = fitlme(lmeAllMice,formulaNumPlaceFields);

formulaMaxPlaceFieldWidth =...
    'maxPlaceFieldWidth~FadPlusorMinus+(1|MouseID)';

lmeMaxPlaceFieldWidth = fitlme(lmeAllMice,formulaMaxPlaceFieldWidth);

plusNumPlaceFields = lmeAllMice{1:176,4};

minusNumPlaceFields = lmeAllMice{177:320,4};

% max place field width
plusMaxPlaceFields = lmeAllMice{1:176,5};

minusMaxPlaceFields = lmeAllMice{177:320,5};