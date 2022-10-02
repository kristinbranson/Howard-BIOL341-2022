function [indicatorLED,isPerStepControl,isRGB] = estimateActivationTiming(expdir,varargin)

[stimulistimingfilestr,protocolfilestr,moviefilestr,offsetRGB,offsetRed,useDataCaptureInfo] = ...
  myparse(varargin,...
  'stimulustimingfilestr','StimulusTimingLog.txt',...
  'protocolfilestr','protocol.mat',...
  'moviefilestr','movie.ufmf',...
  'offsetRGB',[],...
  'offsetRed',[],...
  'useDataCaptureInfo',true);

if isempty(offsetRGB),
  if useDataCaptureInfo,
    offsetRGB = -0.023194;
  else
    offsetRGB = 0.331582;
  end
end
if isempty(offsetRed),
  if useDataCaptureInfo,
    %offsetRed = -0.023194;
    offsetRed = -.123194;
  else
    offsetRed = 0.342450;
  end
end

if useDataCaptureInfo,
  timinginfofile = fullfile(expdir,stimulistimingfilestr);
  [camerastart_datetime,protocolstart_datetime,...
    pulsestart_datetime,~] = ...
    readTimingFile(timinginfofile);
end

protocolfile = fullfile(expdir,protocolfilestr);
if ~exist(protocolfile,'file'),
  % no protocolfile, probably not a chrimson experiment
  indicatorLED = [];
  isPerStepControl = [];
  isRGB = [];
  return;
end
pd = load(protocolfile);
% returns seconds since protocol started
[starttimes,endtimes,intensity,pulsewidths,pulseperiods,...
  stepnums,stepstarts,isRGB] = ...
  parseStimulusProtocol(pd.protocol,'R');

if useDataCaptureInfo,
  isPerStepControl = ~isempty(pulsestart_datetime);

  if isRGB,
    % seconds between camera start and protocol start
    offset = (protocolstart_datetime-camerastart_datetime)*24*3600+offsetRGB;

    starttimes = starttimes + offset;
    endtimes = endtimes + offset;
  else
    assert(numel(pulsestart_datetime) == numel(stepstarts));
    for i = 1:numel(pulsestart_datetime),
      idx = stepnums == i;
      offset = (pulsestart_datetime(i)-camerastart_datetime)*24*3600 - stepstarts(i)+offsetRed;
      starttimes(idx) = starttimes(idx) + offset;
      endtimes(idx) = endtimes(idx) + offset;
    end

  end
else
  isPerStepControl = false;
  if isRGB,
    offset = offsetRGB;
  else
    offset = offsetRed;
  end
  starttimes = starttimes + offset;
  endtimes = endtimes + offset;
end

% which frames do these times correspond to
moviefile = fullfile(expdir,moviefilestr);
headerinfo = ufmf_read_header(moviefile);

startframes = nan(size(starttimes));
endframes = nan(size(starttimes));
fps = (headerinfo.nframes-1)/(headerinfo.timestamps(end)-headerinfo.timestamps(1));
for i = 1:numel(starttimes),
  [dt,fr] = min(abs(headerinfo.timestamps-starttimes(i)));
  if dt > 5*fps,
    fr = nan;
  end
  startframes(i) = fr;
  [dt,fr] = min(abs(headerinfo.timestamps-endtimes(i)));
  if dt > 5*fps,
    fr = nan;
  end
  endframes(i) = fr;
end

% mimic the output of indicator LED detector
indicatordigital = false(1,headerinfo.nframes);
for i = 1:numel(startframes),
  sf = max(1,startframes(i));
  ef = min(headerinfo.nframes,endframes(i));
  indicatordigital(sf:ef-1) = true;
end

indicatorLED = struct;
indicatorLED.StartEndStatus = [isnan(startframes(1)),isnan(endframes(1))];
indicatorLED.startframe = startframes;
indicatorLED.endframe = endframes;
indicatorLED.indicatordigital = indicatordigital;
indicatorLED.starttimes = starttimes;
indicatorLED.endtimes = endtimes;
indicatorLED.intensity = intensity;
indicatorLED.pulsewidths = pulsewidths;
indicatorLED.pulseperiods = pulseperiods;