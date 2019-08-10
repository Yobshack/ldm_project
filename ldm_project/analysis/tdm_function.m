function out = tdm_function(x)

corr_map = corr(x,'rows','pairwise','type','Spearman');

net_corr_light = corr_map(3,1) - corr_map(3,2);
net_corr_dark = corr_map(4,2) - corr_map(4,1);

out = [net_corr_light, net_corr_dark];