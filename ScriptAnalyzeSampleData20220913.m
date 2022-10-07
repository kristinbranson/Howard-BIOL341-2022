%% load in data

load SampleData20221006.mat;
nexps = numel(data.exp);
nflies = nan(1,nexps);
for i = 1:nexps,
  nflies(i) = numel(data.exp(i).fly);
end

ARENARADIUS_MM = 26.689;

exptypes = data.summary.exps.type;
[unique_exptypes,~,exptypeidx] = unique(exptypes);

exptypecolors = lines(numel(exptypes)/2);

minspeed_plot = 0;
maxspeed_plot = 70;

mindist2fly_plot = 0;
maxdist2fly_plot = 15;

fps = median(data.summary.exps.fps);

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
expnum = 1; % which video to plot
flynum = 8; % which fly to plot

hfig = figure(1);
clf;
T = data.exp(expnum).summary.nframes;
plot(data.exp(expnum).fly(flynum).x_mm,data.exp(expnum).fly(flynum).y_mm,'k-');
hold on;
scatter(data.exp(expnum).fly(flynum).x_mm,data.exp(expnum).fly(flynum).y_mm,[],data.exp(expnum).timestamps,'.');
axis equal;
axis(ARENARADIUS_MM*[-1,1,-1,1]);
hcb = colorbar;
hcb.Label.String = 'Time (s)';
xlabel('x-position (mm)');
ylabel('y-position (mm)');
colormap jet;
title(sprintf('Video %d (%s), fly %d (%s)',expnum,data.exp(expnum).summary.type,flynum,data.exp(expnum).fly(flynum).sex),'Interpreter','none');

%% plot the trajectory features for an example fly as a time series
% change these parameters to change which video amd which fly we are plotting

expnum = 1; % which video to plot
flynum = 8; % which fly to plot
prct = 1; % for computing ylims

hfig = figure(2);
clf;
% one set of axes for each trajectory feature
feats = {'x_mm','y_mm','theta_rad','a_mm','b_mm','wing_anglel_rad','wing_angler_rad'};
names = {'x-coord (mm)','y-coord (mm)','orientation (rad)','1/4-major axis length (mm)','1/4-minor axis length (mm)','left wing angle (rad)','right wing angle (rad)'};
nfeat = numel(feats);
hax = gobjects(nfeat,1);
% loop over trajectory features
T = data.exp(expnum).summary.nframes;
for i = 1:nfeat,
  % create the axes
  hax(i) = subplot(nfeat,1,i);
  hold(hax(i),'on');
  feat = data.exp(expnum).fly(flynum).(feats{i});
  fpscurr = data.exp(expnum).summary.fps;
  ylim = [prctile(feat,prct),prctile(feat,100-prct)];
  PlotActivationTimePatch(data.exp(expnum).activation.startframe,data.exp(expnum).activation.endframe,fpscurr,ylim,hax(i));
  % plot the feature for the selected video and fly
  plot((1:T)/fpscurr,feat,'k-');
  xlabel('Time (s)');
  ylabel(names{i},'Interpreter','none');
  set(hax(i),'YLim',ylim);
end
set(hax,'XLim',[1,T]/fps);
% forces the x-axis of all axes to be the same
linkaxes(hax,'x');
title(hax(1),sprintf('Video %d (%s), fly %d (%s)',expnum,data.exp(expnum).summary.type,flynum,data.exp(expnum).fly(flynum).sex),'Interpreter','none');

%% compute the instantaneous speed for each fly

% compute the speed by computing how far the fly moved from one frame to
% the next
for expnum = 1:nexps,
  T = data.exp(expnum).summary.nframes;
  for flynum = 1:numel(data.exp(expnum).fly),
    x = data.exp(expnum).fly(flynum).x_mm;
    y = data.exp(expnum).fly(flynum).y_mm;
    % how far did the fly move between frame t and t + 1 in x?
    dx = x(2:end)-x(1:end-1);
    % how far did the fly move between frame t and t + 1 in y?
    dy = y(2:end)-y(1:end-1);
    % use pythagorean theorem to compute total movement in both x and y
    % normalize by the frames per second to get mm/s
    data.exp(expnum).fly(flynum).speed_mmps = sqrt(dx.^2 + dy.^2) * fps;
  end
end

% statistics over flies
data = ComputeMeanStdErrVideo(data,'speed_mmps');

%% compute the distance to the nearest fly for each fly

% compute the distance between a fly and each other fly, and store the
% minimum distance 
for expnum = 1:nexps,
  T = data.exp(expnum).summary.nframes;
  % loop over main fly
  for mainfly = 1:numel(data.exp(expnum).fly),
    % initialize to nan. uses the trick that the min(something,nan) =
    % something
    data.exp(expnum).fly(mainfly).dist2fly_mm = nan(1,T);
    xmain = data.exp(expnum).fly(mainfly).x_mm;
    ymain = data.exp(expnum).fly(mainfly).y_mm;
    % loop over other flies
    for otherfly = 1:numel(data.exp(expnum).fly),
      % distance to self will always be 0!
      if mainfly == otherfly,
        continue;
      end
      % difference in x-coordinate
      xother = data.exp(expnum).fly(otherfly).x_mm;
      yother = data.exp(expnum).fly(otherfly).y_mm;
      dx = xmain-xother;
      % difference in y-coordinate
      dy = ymain-yother;
      % compute distance using pythagorean theorem
      dcurr = sqrt(dx.^2 + dy.^2);
      % store the minimum between distance to this fly and distance to all
      % previous flies
      data.exp(expnum).fly(mainfly).dist2fly_mm  = min(data.exp(expnum).fly(mainfly).dist2fly_mm,dcurr);
    end
  end
end

% compute stats over all flies
data = ComputeMeanStdErrVideo(data,'dist2fly_mm');

%% plot the flies' speeds for each video

hfig = figure(3);
PlotFeatureOverVideo(data,'speed_mmps',...
  'featlabel','Speed (mm/s)',...
  'plotallflies',false,'plotstderr',true,...
  'maxfeatplot',maxspeed_plot,...
  'minfeatplot',minspeed_plot);

%% plot the distance apart for each fly over the videos

hfig = figure(4);
PlotFeatureOverVideo(data,'dist2fly_mm',...
  'featlabel','Inter-fly dist. (mm)',...
  'plotallflies',false,'plotstderr',true,...
  'maxfeatplot',maxdist2fly_plot,...
  'minfeatplot',mindist2fly_plot);

%% plot the speeds when the LED turns on

hfig = figure(5);
clf;
PlotFeatureOverIntervals(data,'speed_mmps',...
  'featlabel','Speed (mm/s)',...
  'plotallflies',false,'plotstderr',true,...
  'maxfeatplot',maxspeed_plot,...
  'minfeatplot',minspeed_plot,...
  'nframespre_plot',nframespre_speed,...
  'nframespost_plot',nframespost_speed);

%% plot the inter-fly distances when the LED turns on

hfig = figure(6);
clf;
PlotFeatureOverIntervals(data,'dist2fly_mm',...
  'featlabel','Inter-fly dist. (mm)',...
  'plotallflies',false,'plotstderr',true,...
  'maxfeatplot',maxdist2fly_plot,...
  'minfeatplot',mindist2fly_plot,...
  'nframespre_plot',nframespre_dist2fly,...
  'nframespost_plot',nframespost_dist2fly);

%% compute the average speed at LED onset

data = ...
  ComputeMeanStdErrIntervals(data,'speed_mmps','onset',true,...
  'nframespre',nframespre_speed,'nframespost',nframespost_speed);

%% compute the average distance apart at LED onset

data = ...
  ComputeMeanStdErrIntervals(data,'dist2fly_mm','onset',true,...
  'nframespre',nframespre_dist2fly,'nframespost',nframespost_dist2fly);


%% plot the average speed at LED onset

hfig = figure(7);
clf;
hax = PlotMeanFeatureOverIntervals(data,'speed_mmps',...
  'maxfeatplot',maxspeed_plot,...
  'minfeatplot',minspeed_plot);

%set(hax,'yscale','log');

%% plot the average inter-fly distance at LED onset

hfig = figure(8);
clf;
hax = PlotMeanFeatureOverIntervals(data,'dist2fly_mm',...
  'maxfeatplot',maxdist2fly_plot,...
  'minfeatplot',mindist2fly_plot);

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

meanspeedhiston = cell(nexps,1);
meanspeedhistoff = cell(nexps,1);
stdspeedhiston = cell(nexps,1);
stdspeedhistoff = cell(nexps,1);
stderrspeedhiston = cell(nexps,1);
stderrspeedhistoff = cell(nexps,1);
countson = cell(nexps,1);
countsoff = cell(nexps,1);
onintervaltypes = cell(nexps,1);
onintervaltypeidx = cell(nexps,1);

for i = 1:nexps,

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

  speedhistcurr = zeros(numel(bincenters),numel(onstarts),nflies(i));
  isdatacurr = false(numel(onstarts),nflies(i));
  T = size(speed{i},1);

  % histogram each on interval for each fly
  for j = 1:numel(onstarts),
    % which frames to histogram
    t = onstarts(j);
    t0 = min(T,max(1,t+t0_on));
    t1 = min(T,max(1,t+t1_on));

    for k = 1:nflies(i),

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
    ncurr = nnz(idxcurr)*nflies(i);
    speedhistcurr1 = reshape(speedhistcurr(:,idxcurr,:),[numel(bincenters),ncurr]);
    isdatacurr1 = reshape(isdatacurr(idxcurr,:),[1,ncurr]);
    speedhistcurr1 = speedhistcurr1(:,isdatacurr1);
    meanspeedhiston{i}(:,typei) = mean(speedhistcurr1,2);
    stdspeedhiston{i}(:,typei) = std(speedhistcurr1,0,2);
    stderrspeedhiston{i}(:,typei) = stdspeedhiston{i}(:,typei) / sqrt(size(speedhistcurr1,2));
  end

  % histogram each off interval for each fly
  speedhistcurr = zeros(numel(bincenters),numel(offstarts),nflies(i));

  for j = 1:numel(offstarts),
    % which frames to histogram
    t = offstarts(j);
    t0 = min(T,max(1,t+t0_off));
    t1 = min(T,max(1,t+t1_off));

    for k = 1:nflies(i),

      datacurr = speed{i}(t0:t1,k);  
      datacurr = datacurr(~isnan(datacurr));
      isdatacurr(j,k) = numel(datacurr) > (t1-t0)*minfracframeshist;
      speedhistcurr(:,j,k) = histcounts(datacurr,binedges,'Normalization','probability');

    end
  end

  % compute mean, std, stderr
  ncurr = numel(offstarts)*nflies(i);
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

meandist2flyhiston = cell(nexps,1);
meandist2flyhistoff = cell(nexps,1);
stddist2flyhiston = cell(nexps,1);
stddist2flyhistoff = cell(nexps,1);
stderrdist2flyhiston = cell(nexps,1);
stderrdist2flyhistoff = cell(nexps,1);
countson = cell(nexps,1);
countsoff = cell(nexps,1);
onintervaltypes = cell(nexps,1);
onintervaltypeidx = cell(nexps,1);

for i = 1:nexps,

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

  dist2flyhistcurr = zeros(numel(bincenters),numel(onstarts),nflies(i));
  isdatacurr = false(numel(onstarts),nflies(i));

  % histogram each on interval for each fly
  for j = 1:numel(onstarts),
    % which frames to histogram
    t = onstarts(j);
    t0 = min(T-1,max(1,t+t0_on));
    t1 = min(T-1,max(1,t+t1_on));

    for k = 1:nflies(i),

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
    ncurr = nnz(idxcurr)*nflies(i);
    dist2flyhistcurr1 = reshape(dist2flyhistcurr(:,idxcurr,:),[numel(bincenters),ncurr]);
    isdatacurr1 = reshape(isdatacurr(idxcurr,:),[1,ncurr]);
    dist2flyhistcurr1 = dist2flyhistcurr1(:,isdatacurr1);
    meandist2flyhiston{i}(:,typei) = mean(dist2flyhistcurr1,2);
    stddist2flyhiston{i}(:,typei) = std(dist2flyhistcurr1,0,2);
    stderrdist2flyhiston{i}(:,typei) = stddist2flyhiston{i}(:,typei) / sqrt(size(dist2flyhistcurr1,2));
  end

  % histogram each off interval for each fly
  dist2flyhistcurr = zeros(numel(bincenters),numel(offstarts),nflies(i));

  for j = 1:numel(offstarts),
    % which frames to histogram
    t = offstarts(j);
    t0 = min(T-1,max(1,t+t0_off));
    t1 = min(T-1,max(1,t+t1_off));

    for k = 1:nflies(i),

      datacurr = dist2fly{i}(t0:t1,k);  
      datacurr = datacurr(~isnan(datacurr));
      isdatacurr(j,k) = numel(datacurr) > (t1-t0)*minfracframeshist;
      dist2flyhistcurr(:,j,k) = histcounts(datacurr,binedges,'Normalization','probability');

    end
  end

  % compute mean, std, stderr
  ncurr = numel(offstarts)*nflies(i);
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
naxr = ceil(nexps/naxc);
hax = gobjects(naxr,naxc);

% whether to plot all the on intervals, or just the brightest
plotallonintervals = false;

for i = 1:nexps,

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
naxr = ceil(nexps/naxc);
hax = gobjects(naxr,naxc);

% whether to plot all the on intervals, or just the brightest
plotallonintervals = false;

for i = 1:nexps,

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