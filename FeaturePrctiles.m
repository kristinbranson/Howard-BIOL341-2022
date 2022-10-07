function featps = FeaturePrctiles(data,featname,prctiles)

feat = [];
for expnum = 1:numel(data.exp),
  for flynum = 1:numel(data.exp(expnum)),
    idx = ~isnan(data.exp(expnum).fly(flynum).(featname));
    feat = [feat,data.exp(expnum).fly(flynum).(featname)(idx)]; %#ok<AGROW> 
  end
end

featps = prctile(feat,prctiles);
featps(prctiles==100) = max(feat(~isinf(feat)));
featps(prctiles==0) = min(feat(~isinf(feat)));
