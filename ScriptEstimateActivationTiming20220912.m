%% can we figure out the led info without detecting?

for i = 1:numel(expdirs),
  indicatorLED = estimateActivationTiming(expdirs{i});
  id = load(fullfile(allexpdirs{i},'indicatordata.mat'));
  dt = id.indicatorLED.starttimes-indicatorLED.starttimes;
  df = id.indicatorLED.startframe-indicatorLED.startframe;
  [~,expname] = fileparts(expdirs{i});
  fprintf('%s: max diff = %.4f s (%d fr), mean abs diff = %.4f s (%.1f fr), mean diff = %.4f s (%.1f fr)\n',expname,...
    max(abs(dt)),max(abs(df)),mean(abs(dt)),mean(abs(df)),...
    mean(dt),mean(df));
end

%% try to figure out offset

alldts = cell(size(expinfo));
alldfs = cell(size(expinfo));
isPerStepControl = nan(size(expinfo));
isRGB = nan(size(expinfo));
for i = 1:numel(expinfo),
  expdir = expinfo(i).expdir;
  [~,expname] = fileparts(expdir);
  if ~exist(fullfile(expdir,'indicatordata.mat'),'file') || ~exist(fullfile(expdir,'protocol.mat'),'file'),
    continue;
  end
  [indicatorLED,isPerStepControl(i),isRGB(i)] = estimateActivationTiming(expdir,'useDataCaptureInfo',true);
  id = load(fullfile(expdir,'indicatordata.mat'));
  nons_detected = numel(id.indicatorLED.starttimes);
  nons = numel(indicatorLED.starttimes);
  if nons ~= nons_detected,
    fprintf('%d %s: number of on intervals detected (%d) does not match number in protocol (%d)\n',i,expname,nons_detected,nons);
    continue;
  end
  dt = id.indicatorLED.starttimes-indicatorLED.starttimes;
  df = id.indicatorLED.startframe-indicatorLED.startframe;
  alldts{i} = dt;
  alldfs{i} = df;
  fprintf('%d %s: max diff = %.4f s (%d fr), mean abs diff = %.4f s (%.1f fr), mean diff = %.4f s (%.1f fr), perstep = %d, rgb = %d\n',i,expname,...
    max(abs(dt)),max(abs(df)),mean(abs(dt)),mean(abs(df)),...
    mean(dt),mean(df),isPerStepControl(i),isRGB(i));
end

rgbdts = alldts(isRGB == 1 & ~cellfun(@isempty,alldts));
rgboffset = median(cellfun(@mean,rgbdts));
% nons = cellfun(@numel,rgbdts);
% if all(nons==nons(1)),
%   rgbdts = cat(1,rgbdts{:});
%   rgbmeanoffset = mean(rgbdts,1);
% end
reddts = alldts(isRGB == 0 & ~cellfun(@isempty,alldts));
redoffset = median(cellfun(@mean,reddts));
% nons = cellfun(@numel,reddts);
% if all(nons==nons(1)),
%   reddts = cat(1,reddts{:});
%   redmeanoffset = mean(reddts,1);
% end

dfs = alldfs(~cellfun(@isempty,alldfs));
maxerrfr = cellfun(@(x) max(abs(x)),dfs);
meanmaxerrfr = mean(maxerrfr);
maxmaxerrfr = max(maxerrfr);
meanerrfr = mean(abs(cat(2,dfs{:})));
