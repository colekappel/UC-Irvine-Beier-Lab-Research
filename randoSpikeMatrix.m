%vvv start w 5507 vvv this code gets the rando spike matx
function randoMatx = randoSpikeMatrix(eventRatePerBinMatx)

randoMatx=zeros(height(eventRatePerBinMatx),width(eventRatePerBinMatx)); % dimensions of spike matx

for b=1:length(eventRatePerBinMatx)

for a = 1: sum(eventRatePerBinMatx(b,:))
    randoBin = randi(10);
    randoMatx(b,randoBin) = randoMatx(b,randoBin)+1;
end

end
