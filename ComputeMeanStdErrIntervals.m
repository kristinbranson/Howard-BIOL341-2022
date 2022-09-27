function [meanfeatinterval,stderrfeatinterval] = ComputeMeanStdErrIntervals(feat,activation,varargin)

nvideos = numel(feat);

[nframespre,nframespost,tfonset] = ...
  myparse(varargin,...
  'nframespre',30,...
  'nframespost',150,...
  'onset',true);

intervalT = nframespre+nframespost+1;
meanfeatinterval = zeros(intervalT,nvideos);
stderrfeatinterval = zeros(intervalT,nvideos);

if tfonset,
  allts = activation.startframes;
else
  allts = activation.endframes;
end

for i = 1:nvideos,

  T = size(feat{i},1);
  count = zeros(intervalT,1);
  for j = 1:numel(allts{i}),
    t = allts{i}(j);
    t0 = max(1,t-nframespre);
    t1 = min(T,t+nframespost);
    i0 = t0-(t-nframespre) + 1;
    i1 = i0 + (t1-t0); 
    meanfeatinterval(i0:i1,i) = meanfeatinterval(i0:i1,i) + sum(feat{i}(t0:t1,:),2,'omitnan');
    count(i0:i1) = count(i0:i1) + sum(~isnan(feat{i}(t0:t1,:)),2);
  end
  meanfeatinterval(:,i) = meanfeatinterval(:,i) ./ count;

  for j = 1:numel(allts{i}),
    t = allts{i}(j);
    t0 = max(1,t-nframespre);
    t1 = min(T,t+nframespost);
    i0 = t0-(t-nframespre) + 1;
    i1 = i0 + (t1-t0); 
    stderrfeatinterval(i0:i1,i) = stderrfeatinterval(i0:i1,i) + sum((feat{i}(t0:t1,:)-meanfeatinterval(i0:i1,i)).^2,2,'omitnan');
  end
  stderrfeatinterval(:,i) = sqrt(stderrfeatinterval(:,i) ./ count) ./ max(1,sqrt(count));
end