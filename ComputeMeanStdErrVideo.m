function [meanfeat,stderrfeat] = ComputeMeanStdErrVideo(feat)

nvideos = numel(feat);

% average over flies
meanfeat = cell(1,nvideos);
for i = 1:nvideos,
  meanfeat{i} = mean(feat{i},2,'omitnan');
end

% compute the standard deviation over all flies
stderrfeat = cell(1,nvideos);
for i = 1:nvideos,
  ntrajcurr = sum(~isnan(feat{i}),2);
  stderrfeat{i} = std(feat{i},0,2,'omitnan')./sqrt(ntrajcurr);
end