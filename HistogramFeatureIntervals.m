function data = HistogramFeatureIntervals(data,featname,varargin)

nexps = numel(data.exp);

[nframespre,nframespost,tfon,fracfeatname,binsfeatname,deltatname,...
  nbins,binlimits,bincenters] = ...
  myparse(varargin,...
  'nframespre',30,...
  'nframespost',150,...
  'on',true,...
  'fracfeatname','',...
  'binsfeatname','',...
  'deltatname','',...
  'nbins',100,...
  'binlimits',[],...
  'bincenters',[]);

if tfon,
  onname = 'on';
else
  onname = 'off';
end

if isempty(fracfeatname)
  fracfeatname = sprintf('frac_%s_%s',onname,featname);
end
if isempty(binsfeatname)
  binsfeatname = sprintf('bincenters_%s_%s',onname,featname);
end
if isempty(deltatname),
  deltatname = sprintf('deltat_%s_%s',onname,featname);
end

expnums = data.summary.activation.expnum;
if tfon,
  alltstarts = data.summary.activation.startframe;
  alltends = data.summary.activation.endframe;
else
  alltstarts = [1;data.summary.activation.endframe+1];
  alltends = [data.summary.activation.startframe-1;data.summary.exps.nframes];
end

if isempty(bincenters),
  if isempty(binlimits),
    binlimits = FeaturePrctiles(data,featname,[0,100]);
  end
  binedges = linspace(binlimits(1),binlimits(2),nbins+1);
  bincenters = (binedges(1:end-1)+binedges(2:end))/2;
  binedges(1) = -inf; binedges(end) = inf;
else
  binedges = [-inf,(bincenters(1:end-1)+bincenters(2:end))/2,inf];
end


for i = 1:nexps,

  feat = [];
  tstartscurr = alltstarts(expnums==i);
  tendscurr = alltends(expnums==i);

  nframespostcurr = 0;
  for j = 1:numel(tstartscurr),
    tstart = tstartscurr(j);
    tend = tendscurr(j);

    t0 = tstart+nframespre;
    t1 = min(tend,tstart+nframespost);
    if t1 < t0,
      continue;
    end
    nframespostcurr = max(nframespostcurr,t1-tstart);
    for k = 1:numel(data.exp(i).fly),
      featcurr = data.exp(i).fly(k).(featname)(t0:t1);
      feat = [feat,featcurr(~isnan(featcurr))]; %#ok<AGROW> 
    end
  end

  frac = histcounts(feat,binedges,'Normalization','probability');
  data.exp(i).hist.(fracfeatname) = frac;
  data.exp(i).hist.(binsfeatname) = bincenters;
  data.exp(i).hist.(deltatname) = [nframespre,nframespostcurr];

end