function ldm(flyBehavior1,flyBehavior2,flyBehavior3,flyBehavior4,timeBin,numRe)

% Turn biases
x = flyBehavior1;
y = flyBehavior2;

% Turn rates
xx = flyBehavior3;
yy = flyBehavior4;

% mean and variance normalize input vectors
f = @(x) (x - nanmean(x))./nanstd(x);
out.normalization = f;
x = f(x);
y = f(y);

% Null Model = All animals are p = 0.5
xx = xx*timeBin;  % turns in 30 mins;
yy = yy*timeBin;
null = nan(numRe,1);
for jj = 1:numRe
    n1 = binornd(xx,0.5)./xx;
    n1 = (n1 - 0.5)./nanstd(x);    
    n2 = binornd(yy,0.5)./yy;
    n2 = (n2 - 0.5)./nanstd(y);
    null(jj) = nanvar(n1-n2);
end