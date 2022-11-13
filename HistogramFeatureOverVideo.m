function data = HistogramFeatureOverVideo(data,featname,varargin)
% data = HistogramFeatureOverVideo(data,featname,varargin)
% Histograms the feature featname during on/off periods and stores the
% results in data.exp(i).hist.(fracfeatname). 
% Inputs:
% data: struct holding all the data
% featname: name of feature to compute statistics for
% Outputs: 
% data: struct holding all the data, with the new stats added
% fracfeatname: field name for normalized counts for the feature (default:
% ['frac_all_',featname])
% binsfeatname: field name for bin centers (default:
% ['bincenters_all_',featname])
% nbins: number of bins, only used if bincenters is not specified (default: 100)
% binlimits: edges of the first and last bins. if [], minimum and maximum
% of feature are used. only used if bincenters is not specified. (default: [])
% bincenters: centers of bins. if [], uses binlimits and nbins to choose
% bins automatically. (default: [])
nexps = numel(data.exp);

[fracfeatname,binsfeatname,...
  nbins,binlimits,bincenters] = ...
  myparse(varargin,...
  'fracfeatname','',...
  'binsfeatname','',...
  'nbins',100,...
  'binlimits',[],...
  'bincenters',[]);

if isempty(fracfeatname)
  fracfeatname = sprintf('frac_all_%s',featname);
end
if isempty(binsfeatname)
  binsfeatname = sprintf('bincenters_all_%s',featname);
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
  for k = 1:numel(data.exp(i).fly),
    featcurr = data.exp(i).fly(k).(featname);
    feat = [feat,featcurr(~isnan(featcurr))]; %#ok<AGROW>
  end

  frac = histcounts(feat,binedges,'Normalization','probability');
  data.exp(i).hist.(fracfeatname) = frac;
  data.exp(i).hist.(binsfeatname) = bincenters;

end