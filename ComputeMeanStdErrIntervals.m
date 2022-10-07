function data = ComputeMeanStdErrIntervals(data,featname,varargin)

nexps = numel(data.exp);

[nframespre,nframespost,tfonset,meanfeatname,stdfeatname,stderrfeatname,deltatname] = ...
  myparse(varargin,...
  'nframespre',30,...
  'nframespost',150,...
  'onset',true,...
  'meanfeatname','',...
  'stdfeatname','',...
  'stderrfeatname','',...
  'deltatname','');

if tfonset,
  onsetname = 'onset';
else
  onsetname = 'offset';
end

if isempty(meanfeatname)
  meanfeatname = sprintf('mean_%s_%s',onsetname,featname);
end
if isempty(stdfeatname)
  stdfeatname = sprintf('std_%s_%s',onsetname,featname);
end
if isempty(stderrfeatname)
  stderrfeatname = sprintf('stderr_%s_%s',onsetname,featname);
end
if isempty(deltatname),
  deltatname = sprintf('deltat_%s_%s',onsetname,featname);
end

intervalT = nframespre+nframespost+1;

expnums = data.summary.activation.expnum;
if tfonset,
  allts = data.summary.activation.startframe;
else
  allts = activation.endframes;
end

for i = 1:nexps,

  T = numel(data.exp(i).fly(1).(featname));
  tscurr = allts(expnums==i);

  data.exp(i).stat.(meanfeatname) = zeros(1,intervalT);
  data.exp(i).stat.(stdfeatname) = zeros(1,intervalT);
  data.exp(i).stat.(stderrfeatname) = zeros(1,intervalT);
  data.exp(i).stat.(deltatname) = [-nframespre,nframespost];

  count = zeros(1,intervalT);
  for j = 1:numel(tscurr),
    t = tscurr(j);
    t0 = max(1,t-nframespre);
    t1 = min(T,t+nframespost);
    i0 = t0-(t-nframespre) + 1;
    i1 = i0 + (t1-t0);
    for k = 1:numel(data.exp(i).fly),
      feat = data.exp(i).fly(k).(featname)(t0:t1);
      idx = isnan(feat);
      feat(idx) = 0;
      data.exp(i).stat.(meanfeatname)(i0:i1) = data.exp(i).stat.(meanfeatname)(i0:i1) + feat;
      count(i0:i1) = count(i0:i1) + double(~idx);
    end
  end
  data.exp(i).stat.(meanfeatname) = data.exp(i).stat.(meanfeatname) ./ count;

  for j = 1:numel(tscurr),
    t = tscurr(j);
    t0 = max(1,t-nframespre);
    t1 = min(T,t+nframespost);
    i0 = t0-(t-nframespre) + 1;
    i1 = i0 + (t1-t0); 
    feat = data.exp(i).fly(k).(featname)(t0:t1);
    dfeat = (feat-data.exp(i).stat.(meanfeatname)(i0:i1)).^2;
    idx = isnan(dfeat);
    dfeat(idx) = 0;
    data.exp(i).stat.(stdfeatname)(i0:i1) = data.exp(i).stat.(stdfeatname)(i0:i1) + dfeat;
  end
  data.exp(i).stat.(stdfeatname) = sqrt(data.exp(i).stat.(stdfeatname) ./ count);
  data.exp(i).stat.(stdfeatname)(count==0) = nan;
  data.exp(i).stat.(stderrfeatname) = data.exp(i).stat.(stdfeatname) ./ sqrt(count);
  data.exp(i).stat.(stderrfeatname)(count==0) = nan;
end