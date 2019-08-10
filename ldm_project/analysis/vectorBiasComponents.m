function out=vectorBiasComponents(data)

shibireFull=data{1};
controlFull=data{2};

% clear out NaNs
shibireFull(isnan(sum(shibireFull,2)),:)=[];
controlFull(isnan(sum(controlFull,2)),:)=[];

numShib=size(shibireFull,1);
numCtrl=size(controlFull,1);

numBS=5000;
out=zeros(numBS,4);


for i=1:numBS
    
    %resample flies
    shibire=shibireFull(randi(numShib,numShib,1),:);
    control=controlFull(randi(numCtrl,numCtrl,1),:);
    
    % unitize
    shibire=shibire./repmat(sqrt(sum(shibire.^2)),size(shibire,1),1);
    control=control./repmat(sqrt(sum(control.^2)),size(control,1),1);
    
    % compute dot products
    outTemp=[dot(shibire(:,1),shibire(:,3)) dot(shibire(:,2),shibire(:,4)) ...
        dot(control(:,1),control(:,3)) dot(control(:,2),control(:,4))];
    
    outTemp=1-outTemp;
    out(i,:)=outTemp;
    
end

out = [out(:,1) - out(:,3),out(:,2) - out(:,4)];


%plot
plotBool=0;
if plotBool==1
%     figure;
    hold on
    bar(mean(out));
    errorbar(mean(out),std(out),'k.');
    xticks([1 2 3 4]);
    set(gca,'xticklabel',{'shi-L','shi-D','ctr-L','ctr-D'});
end