function out=ldmAlignPair2(a1,a2,i1,i2,v1,v2,metric)
lo_prctile=2;
hi_prctile=98;
a1_5=prctile(a1(:,1),lo_prctile);
a1_95=prctile(a1(:,1),hi_prctile);
a2_5=prctile(a2(:,1),lo_prctile);
a2_95=prctile(a2(:,1),hi_prctile);
a1(:,1)=(a1(:,1)-a1_5)/(a1_95-a1_5);
a2(:,1)=(a2(:,1)-a2_5)/(a2_95-a2_5);
numBins=50;
bins=linspace(-1,2,numBins*3+1);
d1=zeros(length(bins)-1,1);
for i=1:length(bins)-1
    wh=logical(a1(:,1)>=bins(i)) & logical(a1(:,1)<bins(i+1));
    total_intensity=sum(i1(wh));
    total_volume=sum(wh)*v1;
    switch metric
        case 'density'
            d1(i)=total_intensity/total_volume;
        case 'volume'
            d1(i)=total_volume;
        case 'intensity'
            d1(i)=total_intensity;
    end
end
d2=zeros(length(bins)-1,1);
for i=1:length(bins)-1
    wh=logical(a2(:,1)>=bins(i)) & logical(a2(:,1)<bins(i+1));
    total_intensity=sum(i2(wh));
    total_volume=sum(wh)*v2;
    switch metric
        case 'density'
            d2(i)=total_intensity/total_volume;
        case 'volume'
            d2(i)=total_volume;
        case 'intensity'
            d2(i)=total_intensity;
    end
end
kernW=1;
kern=gausswin(kernW);
d1=filter(kern,1,d1);
d2=filter(kern,1,d2);
out.d1=d1;
out.d2=d2;

end