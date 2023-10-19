clear

mouse5507=load('5507Vars.mat');
mouse5508=load('5508Vars.mat');

%box plots for average amplitude
g1 = repmat({'5507 5xFAD+'},length(mouse5507.avAmplitudeForEachCell),1);
g2= repmat({'5508 5xFAD-'},length(mouse5508.avAmplitudeForEachCell),1);
g=[g1;g2];
x=[mouse5507.avAmplitudeForEachCell; mouse5508.avAmplitudeForEachCell];
boxplot(x,g);
[hAmp,pAmp] = ttest2(mouse5507.avAmplitudeForEachCell, mouse5508.avAmplitudeForEachCell);
title('Average Amplitude Comparison')
ylabel('Amplitude (dF/F)')
red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [blue;red];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

%box plots for event rate
figure
gg1 = repmat({'5507 5xFAD+'},length(mouse5507.eventRate),1);
gg2= repmat({'5508 5xFAD-'},length(mouse5508.eventRate),1);
gg=[gg1;gg2];
xx=[mouse5507.eventRate; mouse5508.eventRate];
boxplot(xx,gg);
[hEvent,pEvent] = ttest2(mouse5507.eventRate, mouse5508.eventRate);
title('Event Rate Comparison')
ylabel('Event Rate (Hz)')

red=[255,0,0]/255;
blue =[0,0, 255]/255;
colors = [blue;red];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end