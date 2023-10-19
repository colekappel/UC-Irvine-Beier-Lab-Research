function score = infoScore(eventRatePerBinMatx, Pi)
%Pi is probability mouse in ith bin (time mouse in ith bin/ whole time, 10
%elements in Pi vector)
%lambda is average event rate for each cell across trial

score=[];
for a =1:height(eventRatePerBinMatx)
    theInfoScore=0;
    lambda = mean(eventRatePerBinMatx(a,:)); %average event rate per cell across trial (lambda in formula)
    for b=1:10
        if eventRatePerBinMatx(a,b) ==0
            theInfoScore=theInfoScore+0;
        else
            theInfoScore=theInfoScore+ Pi(b)*(eventRatePerBinMatx(a,b)/lambda)*log2(eventRatePerBinMatx(a,b)/lambda);
        end
    end
    score(a)=theInfoScore;
end

end