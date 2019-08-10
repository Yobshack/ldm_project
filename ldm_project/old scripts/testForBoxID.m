

%cell2array(allDataCell{1,9})

x = allDataCell{1,9}


x = cellfun(@(x) str2double(x),allDataCell(:,9),'un',false);
y = cat(1,x{:});


boxes = cell2mat(allDataCell(:,7));

[a,~,c] = unique(y(:,1:4),'rows');

for ii = 1:size(a,1)
bxs = boxes(c == ii);
mins = y(c == ii,5);

notNans = ~isnan(bxs);
bxs = bxs(notNans);
mins = mins(notNans);

[d,~,f] = unique(bxs);
d
hold on
dat = grpstats(mins,f,'mean');
plot(1:numel(d),dat,'Marker','o','MarkerFaceColor','black')
text(1:numel(d),dat,num2str(a(ii,:)))
 
end

shg

%%


dateNums = cellfun(@(x) datenum(x),adca(:,9));
boxes = cell2mat(adca(:,7));
plot(dateNums,boxes,'Marker','o','LineStyle','none')
datetick('x')
shg
