%% load in data

load SampleData20220915.mat;
nvideos = numel(data);
nfeat = size(data{1},2);
ntrajs = nan(1,nvideos);
for i = 1:nvideos,
  ntrajs(i) = size(data{i},3);
end

XIDX = 1;
YIDX = 2;
ORIIDX = 3;
WINGANGLEL_IDX = 4;
WINGANGLER_IDX = 5;

ARENARADIUS_MM = 26.689;

genotypes = regexp(expnames,'^(.*)_Rig','tokens','once');
genotypes = [genotypes{:}];
[unique_genotypes,~,genotypeidx] = unique(genotypes);

genotypecolors = lines(numel(genotypes)/2);

minspeed_plot = 0;
maxspeed_plot = 70;

mindist2fly_plot = 0;
maxdist2fly_plot = 15;

nframespre_speed = 1*fps;
nframespost_speed = 5*fps;

nframespre_dist2fly = 5*fps;
nframespost_dist2fly = 30*fps;


% data: tracked positions of the flies in sample videos.
% data is a cell with an entry for each video. data{i} is a matrix of size 
% T x 5 x ntraj. ntraj is the number of trajectories. This ideally would be
% the number of flies in the video, but might be more if the tracker loses
% track of a fly for a bit. In this case, you will see a trajectory that
% has data for the first t frames, and is nan for the rest, then a second
% trajectory that has data only for the last T - t + 1 frames. 
% data{i}(:,1,:) is the x-coordinate of the centroid of the fly in
% millimeters. 
% data{i}(:,2,:) is the y-coordinate of the centroid of the fly in
% millimeters. 
% data{i}(:,3,:) is the orientation of the fly in radians.
% data{i}(:,4,:) is the angle of the fly's left wing in radians.
% data{i}(:,5,:) is the angle of the fly's right wing in radians.
% 
% names: 1 x 5 cell with names of the tracking features in data. 
%
% sex: whether each fly is female ('f') or male ('m'). sex is a cell with
% an entry for each video. sex{i}(j) is the sex of fly j for video i. 
% 
% activation: information about when the red LED panel state during the
% experiment. struct with the following fields:
% activation.startframes: when the LEDs turn on in frames.
% cell with an entry for each video. 
% activation.endframes: when the LEDs turn off in frames. cell with an
% entry for each video.
% activation.intensities: how bright the LEDs are during the stimulus
% period. cell with an entry for each video. 
% activation.pulsewidths, activation.pulseperiods: in some cases, the LEDs
% are strobed with an on time of pulsewidth per pulseperiod, measured in
% milliseconds. If the LED is constantly on, the pulsewidth and pulseperiod
% will be the same. cell with an entry for each video.
% 
% expnames: cell with an entry for each video containing metadata about the
% video, including fly genotype, when the data was collected, which rig it
% was collected on. 
% 
% fps: frames per second. scalar. 

%% plot an example fly's centroid trajectory

% change these parameters to change which video, which fly, and which
% frames we are plotting
videoi = 1; % which video to plot
fly = 8; % which fly to plot

hfig = figure(1);
clf;
T = size(data{videoi},1);
plot(data{videoi}(:,XIDX,fly),data{videoi}(:,YIDX,fly),'k-');
hold on;
scatter(data{videoi}(:,XIDX,fly),data{videoi}(:,YIDX,fly),[],(1:T)/fps,'.');
axis equal;
axis(ARENARADIUS_MM*[-1,1,-1,1]);
hcb = colorbar;
hcb.Label.String = 'Time (s)';
xlabel('x-position (mm)');
ylabel('y-position (mm)');
colormap jet;
title(sprintf('Video %d (%s), fly %d (%s)',videoi,expnames{videoi},fly,sex{videoi}(fly)),'Interpreter','none');

%% plot the trajectory features for an example fly as a time series
% change these parameters to change which video amd which fly we are plotting

videoi = 1; % which video to plot
fly = 2; % which fly to plot

hfig = figure(2);
clf;
% one set of axes for each trajectory feature
hax = gobjects(nfeat);
% loop over trajectory features
for i = 1:nfeat,
  T = size(data{videoi},1);
  % create the axes
  hax(i) = subplot(nfeat,1,i);
  hold(hax(i),'on');
  ylim = [min(data{videoi}(:,i,fly)),max(data{videoi}(:,i,fly))];
  PlotActivationTimePatch(activation.startframes{videoi},activation.endframes{videoi},fps,ylim,hax(i));
  % plot the feature for the selected video and fly
  plot((1:T)/fps,data{videoi}(:,i,fly),'k-');
  xlabel('Time (s)');
  ylabel(names{i},'Interpreter','none');
end
set(hax,'XLim',[t0,t1]/fps);
% forces the x-axis of all axes to be the same
linkaxes(hax,'x');
title(hax(1),sprintf('Video %d (%s), fly %d (%s)',videoi,expnames{videoi},fly,sex{videoi}(fly)),'Interpreter','none');

%% compute the instantaneous speed for each fly

% compute the speed by computing how far the fly moved from one frame to
% the next
speed = cell(1,nvideos);
for i = 1:nvideos,
  T = size(data{i},1);
  % how far did the fly move between frame t and t + 1 in x?
  dx = data{i}(2:end,XIDX,:)-data{i}(1:end-1,XIDX,:);
  % how far did the fly move between frame t and t + 1 in y?
  dy = data{i}(2:end,YIDX,:)-data{i}(1:end-1,YIDX,:);
  % use pythagorean theorem to compute total movement in both x and y
  % normalize by the frames per second to get mm/s
  speed{i} = sqrt(dx.^2 + dy.^2) * fps;
  % reshape to T-1 x ntrajs(i)
  speed{i} = reshape(speed{i},[T-1,ntrajs(i)]);
end

% statistics over flies
[meanspeed,stderrspeed] = ComputeMeanStdErrVideo(speed);

%% compute the distance to the nearest fly

% compute the distance between a fly and each other fly, and store the
% minimum distance 
dist2fly = cell(1,nvideos);
for i = 1:nvideos,
  T = size(data{i},1);
  % initialize to nan. uses the trick that the min(something,nan) =
  % something
  dist2fly{i} = nan(T,ntrajs(i));
  % loop over main fly
  for fly1 = 1:ntrajs(i),
    % loop over other flies
    for fly2 = 1:ntrajs(i),
      % distance to self will always be 0!
      if fly1 == fly2,
        continue;
      end
      % difference in x-coordinate
      dx = data{i}(:,XIDX,fly1)-data{i}(:,XIDX,fly2);
      % difference in y-coordinate
      dy = data{i}(:,YIDX,fly1)-data{i}(:,YIDX,fly2);
      % compute distance using pythagorean theorem
      dcurr = sqrt(dx.^2 + dy.^2);
      % store the minimum between distance to this fly and distance to all
      % previous flies
      dist2fly{i}(:,fly1) = min(dist2fly{i}(:,fly1),dcurr);
    end
  end
end

% compute stats over all flies
[meandist2fly,stderrdist2fly] = ComputeMeanStdErrVideo(dist2fly);

%% plot the flies' speeds for each video

hfig = figure(3);
PlotFeatureOverVideo(speed,meanspeed,activation,fps,...
  'featlabel','Speed (mm/s)',...
  'plotallflies',false,'plotstderr',true,...
  'genotypeidx',genotypeidx,...
  'stderrfeat',stderrspeed,...
  'expnames',expnames,...
  'maxfeatplot',maxspeed_plot,...
  'minfeatplot',minspeed_plot);

%% plot the distance apart for each fly over the videos

hfig = figure(4);
PlotFeatureOverVideo(dist2fly,meandist2fly,activation,fps,...
  'featlabel','Inter-fly dist. (mm)',...
  'plotallflies',false,'plotstderr',true,...
  'genotypeidx',genotypeidx,...
  'stderrfeat',stderrdist2fly,...
  'expnames',expnames,...
  'maxfeatplot',maxdist2fly_plot,...
  'minfeatplot',mindist2fly_plot);

%% plot the speeds when the LED turns on

hfig = figure(5);
clf;
PlotFeatureOverIntervals(speed,meanspeed,activation,fps,...
  'featlabel','Speed (mm/s)',...
  'plotallflies',false,'plotstderr',true,...
  'genotypeidx',genotypeidx,...
  'stderrfeat',stderrspeed,...
  'maxfeatplot',maxspeed_plot,...
  'minfeatplot',minspeed_plot,...
  'nframespre_plot',nframespre_speed,...
  'nframespost_plot',nframespost_speed,...
  'expnames',expnames);

%% plot the inter-fly distances when the LED turns on

hfig = figure(6);
clf;
PlotFeatureOverIntervals(dist2fly,meandist2fly,activation,fps,...
  'featlabel','Inter-fly dist. (mm)',...
  'plotallflies',false,'plotstderr',true,...
  'genotypeidx',genotypeidx,...
  'stderrfeat',stderrdist2fly,...
  'maxfeatplot',maxdist2fly_plot,...
  'minfeatplot',mindist2fly_plot,...
  'nframespre_plot',nframespre_dist2fly,...
  'nframespost_plot',nframespost_dist2fly,...
  'expnames',expnames);

%% compute the average speed at LED onset

[meanspeed_interval,stderrspeed_interval] = ...
  ComputeMeanStdErrIntervals(speed,activation,...
  'nframespre',nframespre_speed,'nframespost',nframespost_speed);

%% compute the average distance apart at LED onset

[meandist2fly_interval,stderrdist2fly_interval] = ...
  ComputeMeanStdErrIntervals(dist2fly,activation,...
  'nframespre',nframespre_dist2fly,'nframespost',nframespost_dist2fly);


%% plot the average speed at LED onset

hfig = figure(7);
clf;
hax = PlotMeanFeatureOverIntervals(meanspeed_interval,...
  stderrspeed_interval,fps,...
  nframespre_speed,nframespost_speed,...
  'featlabel','Speed (mm/s)',...
  'genotypeidx',genotypeidx,...
  'maxfeatplot',maxspeed_plot,...
  'minfeatplot',minspeed_plot,...
  'expnames',expnames);

%set(hax,'yscale','log');

%% plot the average inter-fly distance at LED onset

hfig = figure(8);
clf;
hax = PlotMeanFeatureOverIntervals(meandist2fly_interval,...
  stderrdist2fly_interval,fps,...
  nframespre_dist2fly,nframespost_dist2fly,...
  'featlabel','Inter-fly dist. (mm)',...
  'genotypeidx',genotypeidx,...
  'maxfeatplot',maxdist2fly_plot,...
  'minfeatplot',mindist2fly_plot,...
  'expnames',expnames);

%% histogram speed during on and off periods

t0_on = 0*fps; % how long after lights on to include
t1_on = 10*fps; % maximum time after lights on to include. set to inf to include to the end of the period. 
t0_off = 20*fps; % how long after lights off to include
t1_off = inf; % maximum time after lights off to include. set to inf to include to the end of the period. 
minfracframeshist = .5; % ignore trajectories that have too few non-nan frames

maxspeedplot = 100; % histogram bin limits -- values outside will be placed in the end bin
nbins = 20; % number of bins

% histogram bin locations
binedges = logspace(-1,log10(maxspeedplot),nbins+1);
bincenters = (binedges(1:end-1)+binedges(2:end))/2;
% plot the fraction > maxspeedplot in a last bin
bw = binedges(2)/binedges(1);
bincenters(end+1) = bincenters(end)*bw;
binedges(end+1) = inf;
binedges(1) = 0;

meanspeedhiston = cell(nvideos,1);
meanspeedhistoff = cell(nvideos,1);
stdspeedhiston = cell(nvideos,1);
stdspeedhistoff = cell(nvideos,1);
stderrspeedhiston = cell(nvideos,1);
stderrspeedhistoff = cell(nvideos,1);
countson = cell(nvideos,1);
countsoff = cell(nvideos,1);
onintervaltypes = cell(nvideos,1);
onintervaltypeidx = cell(nvideos,1);

for i = 1:nvideos,

  % when do on and off periods start
  onstarts = activation.startframes{i};
  onends = activation.endframes{i};
  offstarts = [1,activation.endframes{i}+1];
  offends = [activation.startframes{i}-1,T];

  % group together data from stimuli that are the same
  [onintervaltypes{i},~,onintervaltypeidx{i}] = unique([activation.intensities{i}(:),activation.pulsewidths{i}(:),activation.pulseperiods{i}(:)],'rows');
  non = size(onintervaltypes{i},1);
  % all off stimuli are the same
  noff = 1;

  speedhistcurr = zeros(numel(bincenters),numel(onstarts),ntrajs(i));
  isdatacurr = false(numel(onstarts),ntrajs(i));
  T = size(speed{i},1);

  % histogram each on interval for each fly
  for j = 1:numel(onstarts),
    % which frames to histogram
    t = onstarts(j);
    t0 = min(T,max(1,t+t0_on));
    t1 = min(T,max(1,t+t1_on));

    for k = 1:ntrajs(i),

      datacurr = speed{i}(t0:t1,k);  
      datacurr = datacurr(~isnan(datacurr));
      isdatacurr(j,k) = numel(datacurr) > (t1-t0)*minfracframeshist;
      speedhistcurr(:,j,k) = histcounts(datacurr,binedges,'Normalization','probability');

    end
  end

  meanspeedhiston{i} = zeros(numel(bincenters),non);
  stdspeedhiston{i} = zeros(numel(bincenters),non);
  stderrspeedhiston{i} = zeros(numel(bincenters),non);

  % compute mean, std, stderr
  for typei = 1:non,
    idxcurr = onintervaltypeidx{i} == typei;
    ncurr = nnz(idxcurr)*ntrajs(i);
    speedhistcurr1 = reshape(speedhistcurr(:,idxcurr,:),[numel(bincenters),ncurr]);
    isdatacurr1 = reshape(isdatacurr(idxcurr,:),[1,ncurr]);
    speedhistcurr1 = speedhistcurr1(:,isdatacurr1);
    meanspeedhiston{i}(:,typei) = mean(speedhistcurr1,2);
    stdspeedhiston{i}(:,typei) = std(speedhistcurr1,0,2);
    stderrspeedhiston{i}(:,typei) = stdspeedhiston{i}(:,typei) / sqrt(size(speedhistcurr1,2));
  end

  % histogram each off interval for each fly
  speedhistcurr = zeros(numel(bincenters),numel(offstarts),ntrajs(i));

  for j = 1:numel(offstarts),
    % which frames to histogram
    t = offstarts(j);
    t0 = min(T,max(1,t+t0_off));
    t1 = min(T,max(1,t+t1_off));

    for k = 1:ntrajs(i),

      datacurr = speed{i}(t0:t1,k);  
      datacurr = datacurr(~isnan(datacurr));
      isdatacurr(j,k) = numel(datacurr) > (t1-t0)*minfracframeshist;
      speedhistcurr(:,j,k) = histcounts(datacurr,binedges,'Normalization','probability');

    end
  end

  % compute mean, std, stderr
  ncurr = numel(offstarts)*ntrajs(i);
  speedhistcurr1 = reshape(speedhistcurr,[numel(bincenters),ncurr]);
  isdatacurr1 = reshape(isdatacurr,[1,ncurr]);
  speedhistcurr1 = speedhistcurr1(:,isdatacurr1);
  meanspeedhistoff{i} = mean(speedhistcurr1,2);
  stdspeedhistoff{i} = std(speedhistcurr1,0,2);
  stderrspeedhistoff{i} = stdspeedhistoff{i} / sqrt(size(speedhistcurr1,2));

end

%% histogram dist2fly during on and off periods

t0_on = 0*fps; % how long after lights on to include
t1_on = 30*fps; % maximum time after lights on to include. set to inf to include to the end of the period. 
t0_off = 20*fps; % how long after lights off to include
t1_off = inf; % maximum time after lights off to include. set to inf to include to the end of the period. 
minfracframeshist = .5; % ignore trajectories that have too few non-nan frames

mindistplot = 2;
maxdistplot = 2*ARENARADIUS_MM; % histogram bin limits -- values outside will be placed in the end bin
nbins = 20; % number of bins

% histogram bin locations
binedges = logspace(log10(mindistplot),log10(maxdistplot),nbins+1);
bincenters = (binedges(1:end-1)+binedges(2:end))/2;
% plot the fraction > maxdistplot in a last bin
bw = binedges(2)/binedges(1);
bincenters(end+1) = bincenters(end)*bw;
binedges(end+1) = inf;
binedges(1) = 0;

meandist2flyhiston = cell(nvideos,1);
meandist2flyhistoff = cell(nvideos,1);
stddist2flyhiston = cell(nvideos,1);
stddist2flyhistoff = cell(nvideos,1);
stderrdist2flyhiston = cell(nvideos,1);
stderrdist2flyhistoff = cell(nvideos,1);
countson = cell(nvideos,1);
countsoff = cell(nvideos,1);
onintervaltypes = cell(nvideos,1);
onintervaltypeidx = cell(nvideos,1);

for i = 1:nvideos,

  % when do on and off periods start
  onstarts = activation.startframes{i};
  onends = activation.endframes{i};
  offstarts = [1,activation.endframes{i}+1];
  offends = [activation.startframes{i}-1,T];

  % group together data from stimuli that are the same
  [onintervaltypes{i},~,onintervaltypeidx{i}] = unique([activation.intensities{i}(:),activation.pulsewidths{i}(:),activation.pulseperiods{i}(:)],'rows');
  non = size(onintervaltypes{i},1);
  % all off stimuli are the same
  noff = 1;

  dist2flyhistcurr = zeros(numel(bincenters),numel(onstarts),ntrajs(i));
  isdatacurr = false(numel(onstarts),ntrajs(i));

  % histogram each on interval for each fly
  for j = 1:numel(onstarts),
    % which frames to histogram
    t = onstarts(j);
    t0 = min(T-1,max(1,t+t0_on));
    t1 = min(T-1,max(1,t+t1_on));

    for k = 1:ntrajs(i),

      datacurr = dist2fly{i}(t0:t1,k);  
      datacurr = datacurr(~isnan(datacurr));
      isdatacurr(j,k) = numel(datacurr) > (t1-t0)*minfracframeshist;
      dist2flyhistcurr(:,j,k) = histcounts(datacurr,binedges,'Normalization','probability');

    end
  end

  meandist2flyhiston{i} = zeros(numel(bincenters),non);
  stddist2flyhiston{i} = zeros(numel(bincenters),non);
  stderrdist2flyhiston{i} = zeros(numel(bincenters),non);

  % compute mean, std, stderr
  for typei = 1:non,
    idxcurr = onintervaltypeidx{i} == typei;
    ncurr = nnz(idxcurr)*ntrajs(i);
    dist2flyhistcurr1 = reshape(dist2flyhistcurr(:,idxcurr,:),[numel(bincenters),ncurr]);
    isdatacurr1 = reshape(isdatacurr(idxcurr,:),[1,ncurr]);
    dist2flyhistcurr1 = dist2flyhistcurr1(:,isdatacurr1);
    meandist2flyhiston{i}(:,typei) = mean(dist2flyhistcurr1,2);
    stddist2flyhiston{i}(:,typei) = std(dist2flyhistcurr1,0,2);
    stderrdist2flyhiston{i}(:,typei) = stddist2flyhiston{i}(:,typei) / sqrt(size(dist2flyhistcurr1,2));
  end

  % histogram each off interval for each fly
  dist2flyhistcurr = zeros(numel(bincenters),numel(offstarts),ntrajs(i));

  for j = 1:numel(offstarts),
    % which frames to histogram
    t = offstarts(j);
    t0 = min(T-1,max(1,t+t0_off));
    t1 = min(T-1,max(1,t+t1_off));

    for k = 1:ntrajs(i),

      datacurr = dist2fly{i}(t0:t1,k);  
      datacurr = datacurr(~isnan(datacurr));
      isdatacurr(j,k) = numel(datacurr) > (t1-t0)*minfracframeshist;
      dist2flyhistcurr(:,j,k) = histcounts(datacurr,binedges,'Normalization','probability');

    end
  end

  % compute mean, std, stderr
  ncurr = numel(offstarts)*ntrajs(i);
  dist2flyhistcurr1 = reshape(dist2flyhistcurr,[numel(bincenters),ncurr]);
  isdatacurr1 = reshape(isdatacurr,[1,ncurr]);
  dist2flyhistcurr1 = dist2flyhistcurr1(:,isdatacurr1);
  meandist2flyhistoff{i} = mean(dist2flyhistcurr1,2);
  stddist2flyhistoff{i} = std(dist2flyhistcurr1,0,2);
  stderrdist2flyhistoff{i} = stddist2flyhistoff{i} / sqrt(size(dist2flyhistcurr1,2));

end

%% plot speed histograms

hfig = figure(9);
clf;

% create set of axes for each video
naxc = 2;
naxr = ceil(nvideos/naxc);
hax = gobjects(naxr,naxc);

% whether to plot all the on intervals, or just the brightest
plotallonintervals = false;

for i = 1:nvideos,

  % choose the brightest on interval to plot
  if plotallonintervals,
    idxonplot = 1:size(meanspeedhiston{i},2);
  else
    [~,idxonplot] = max(onintervaltypes{i}(:,1));
  end
  idxoffplot = 1;

  non = numel(idxonplot);
  noff = numel(idxoffplot);
  if plotallonintervals,
    colorson = lines(non);
  else
    colorson = [0.6350    0.0780    0.1840];
  end
  colorsoff = zeros(noff,3);

  hax(i) = subplot(naxr,naxc,i);
  hold(hax(i),'on');

  h = gobjects(non+noff,1);
  legs = cell(non+noff,1);
  for j = 1:non,
    plot(bincenters,meanspeedhiston{i}(:,idxonplot(j))-stderrspeedhiston{i}(:,idxonplot(j)),'-','Color',colorson(j,:));
    plot(bincenters,meanspeedhiston{i}(:,idxonplot(j))+stderrspeedhiston{i}(:,idxonplot(j)),'-','Color',colorson(j,:));
    h(j) = plot(bincenters,meanspeedhiston{i}(:,idxonplot(j)),'-','Color',colorson(j,:),'LineWidth',2);
    if plotallonintervals,
      legs{j} = sprintf('On %d, %d/%d',onintervaltypes{i}(idxonplot(j),:));
    else
      legs{j} = 'On';
    end
  end
  for j = 1:noff,
    plot(bincenters,meanspeedhistoff{i}(:,idxoffplot(j))-stderrspeedhistoff{i}(:,idxoffplot(j)),'-','Color',colorsoff(j,:));
    plot(bincenters,meanspeedhistoff{i}(:,idxoffplot(j))+stderrspeedhistoff{i}(:,idxoffplot(j)),'-','Color',colorsoff(j,:));

    h(j+non) = plot(bincenters,meanspeedhistoff{i}(:,idxoffplot(j)),'-','Color',colorsoff(j,:),'LineWidth',2);
    legs{j+non} = 'Off';
  end
  if plotallonintervals,
    legend(h,legs,'Location','EastOutside');
  elseif i == 1,
    legend(h,legs);
  end
  title(sprintf('Video %d: %s',i,expnames{i}),'Interpreter','none');

end
xlabel(hax(end),'Speed (mm/s)');
ylabel(hax(end),'Fraction of frames');
set(hax,'XLim',bincenters([1,end]),'xscale','log');
linkaxes(hax);


%% plot dist2fly histograms

hfig = figure(10);
clf;

% create set of axes for each video
naxc = 2;
naxr = ceil(nvideos/naxc);
hax = gobjects(naxr,naxc);

% whether to plot all the on intervals, or just the brightest
plotallonintervals = false;

for i = 1:nvideos,

  % choose the brightest on interval to plot
  if plotallonintervals,
    idxonplot = 1:size(meandist2flyhiston{i},2);
  else
    [~,idxonplot] = max(onintervaltypes{i}(:,1));
  end
  idxoffplot = 1;

  non = numel(idxonplot);
  noff = numel(idxoffplot);
  if plotallonintervals,
    colorson = lines(non);
  else
    colorson = [0.6350    0.0780    0.1840];
  end
  colorsoff = zeros(noff,3);

  hax(i) = subplot(naxr,naxc,i);
  hold(hax(i),'on');

  h = gobjects(non+noff,1);
  legs = cell(non+noff,1);
  for j = 1:non,
    plot(bincenters,meandist2flyhiston{i}(:,idxonplot(j))-stderrdist2flyhiston{i}(:,idxonplot(j)),'-','Color',colorson(j,:));
    plot(bincenters,meandist2flyhiston{i}(:,idxonplot(j))+stderrdist2flyhiston{i}(:,idxonplot(j)),'-','Color',colorson(j,:));
    h(j) = plot(bincenters,meandist2flyhiston{i}(:,idxonplot(j)),'-','Color',colorson(j,:),'LineWidth',2);
    if plotallonintervals,
      legs{j} = sprintf('On %d, %d/%d',onintervaltypes{i}(idxonplot(j),:));
    else
      legs{j} = 'On';
    end
  end
  for j = 1:noff,
    plot(bincenters,meandist2flyhistoff{i}(:,idxoffplot(j))-stderrdist2flyhistoff{i}(:,idxoffplot(j)),'-','Color',colorsoff(j,:));
    plot(bincenters,meandist2flyhistoff{i}(:,idxoffplot(j))+stderrdist2flyhistoff{i}(:,idxoffplot(j)),'-','Color',colorsoff(j,:));

    h(j+non) = plot(bincenters,meandist2flyhistoff{i}(:,idxoffplot(j)),'-','Color',colorsoff(j,:),'LineWidth',2);
    legs{j+non} = 'Off';
  end
  if plotallonintervals,
    legend(h,legs,'Location','EastOutside');
  elseif i == 1,
    legend(h,legs);
  end
  title(sprintf('Video %d: %s',i,expnames{i}),'Interpreter','none');

end
xlabel(hax(end),'Inter-fly distance (mm)');
ylabel(hax(end),'Fraction of frames');
set(hax,'XLim',bincenters([1,end]),'xscale','log');
linkaxes(hax);