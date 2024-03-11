function sparsity = sparsityFunc(eventRatePerBinMatx,Pi)
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
end