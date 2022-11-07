function data = HistogramFeatureIntervals(data,featname,varargin)
% data = HistogramFeatureIntervals(data,featname,varargin)
% Histograms the feature featname during on/off periods and stores the
% results in data.exp(i).hist.(fracfeatname). 
% Inputs:
% data: struct holding all the data
% featname: name of feature to compute statistics for
% Outputs: 
% data: struct holding all the data, with the new stats added
% Optional inputs:
% nframespre: how much after the lights on/off to skip (can be negative) (default: 30)
% nframespost: how much before the end of the lights on/off period to skip (default: 30)
% on: whether to compute interval around lights on (true) or lights off
% (false) (default: true). 
% For the following names, if computing offset periods, replace 'on' with
% 'off'.
% fracfeatname: field name for normalized counts for the feature (default:
% ['frac_on_',featname])
% binsfeatname: field name for bin centers (default:
% ['bincenters_on_',featname])
% deltatname: field name to store interval information in (default:
% ['deltat_on_',featname]
% nbins: number of bins, only used if bincenters is not specified (default: 100)
% binlimits: edges of the first and last bins. if [], minimum and maximum
% of feature are used. only used if bincenters is not specified. (default: [])
% bincenters: centers of bins. if [], uses binlimits and nbins to choose
% bins automatically. (default: [])
nexps = numel(data.exp);

[nframespre,nframespost,tfon,fracfeatname,binsfeatname,deltatname,...
  nbins,binlimits,bincenters] = ...
  myparse(varargin,...
  'nframespre',30,...
  'nframespost',30,...
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