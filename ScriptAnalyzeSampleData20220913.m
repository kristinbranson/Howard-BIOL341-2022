%% Script to compare behavior of flies

usesampledata = false;

%% set up paths for the expriment you want to look at and load it in

if usesampledata,

  sampledatafile = 'C:\Code\Howard-BIOL341-2022\SampleData20221006.mat';
  load(sampledatafile);

else

  % where the data directories are (change)
  % rootdatadir = '/groups/branson/bransonlab/alice/temp_howard';
  % rootdatadir = 'E:\BIOL341\GoogleDrive';
  rootdatadir = 'C:\Code\Howard-BIOL341-2022\sample_processed_data\';

  % Names of single experiments by directory name:
  expnames = {
    'CsChrSocial3_aIPg_RigF_20220923T163615'
    'HU_Back_Ctrl_RigE_20220927T172346'
    'HU_Back_Ctrl_RigF_20220927T172138'
    };

  % combine them to make a full path
  expdirs = cell(size(expnames));
  for i = 1:numel(expnames)
    expdirs{i} = fullfile(rootdatadir,expnames{i});
    if ~exist(expdirs{i},'dir'),
      error('Directory %s does not exist',expdirs{i});
    end
  end
  fprintf('Here is the first full exp dir:\n%s\n', expdirs{i});

  data = LoadTracking(expdirs);

end

nexps = numel(data.exp);
nflies = nan(1,nexps);
for i = 1:nexps,
  nflies(i) = numel(data.exp(i).fly);
end

%% set some parameters

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


%% plot an example fly's centroid trajectory

% change these parameters to change which video, which fly, and which
% frames we are plotting
expnum = 1; % which video to plot
flynum = 7; % which fly to plot

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
flynum = 7; % which fly to plot
prct = 1; % for computing ylims

hfig = figure(2);
clf;
% one set of axes for each trajectory feature
feats = {'x_mm','y_mm','theta_rad','a_mm','b_mm','wing_anglel_rad','wing_angler_rad'};
names = {'x-coord (mm)','y-coord (mm)','orientation (rad)','1/4-major axis (mm)','1/4-minor axis (mm)','left wing angle (rad)','right wing angle (rad)'};
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
  'minfeatplot',minspeed_plot,...
  'featlabel','Speed (mm/s)');


%set(hax,'yscale','log');

%% plot the average inter-fly distance at LED onset

hfig = figure(8);
clf;
hax = PlotMeanFeatureOverIntervals(data,'dist2fly_mm',...
  'maxfeatplot',maxdist2fly_plot,...
  'minfeatplot',mindist2fly_plot,...
  'featlabel','Inter-fly dist. (mm)');


%% histogram speed during on and off periods

t0_on = round(.25*fps); % how long after lights on to include
t1_on = inf; % maximum time after lights on to include. set to inf to include to the end of the period. 
t0_off = round(.25*fps); % how long after lights off to include
t1_off = inf; % maximum time after lights off to include. set to inf to include to the end of the period. 

nbins = 25;
binlimits_speed = [0,maxspeed_plot];

data = HistogramFeatureIntervals(data,'speed_mmps',...
  'nframespre',t0_on,...
  'nframespost',t1_on,...
  'on',true,...
  'nbins',nbins,...
  'binlimits',binlimits_speed);

data = HistogramFeatureIntervals(data,'speed_mmps',...
  'nframespre',t0_off,...
  'nframespost',t1_off,...
  'on',false,...
  'nbins',nbins,...
  'binlimits',binlimits_speed);

hfig = figure(9);
clf;
hax = PlotHistogramOnVsOff(data,'speed_mmps',...
  'featlabel','Speed (mm/s)');


%% histogram dist2fly during on and off periods

t0_on = round(.25*fps); % how long after lights on to include
t1_on = inf; % maximum time after lights on to include. set to inf to include to the end of the period. 
t0_off = round(.25*fps); % how long after lights off to include
t1_off = inf; % maximum time after lights off to include. set to inf to include to the end of the period. 

nbins = 25;
binlimits_speed = [0,maxspeed_plot];

data = HistogramFeatureIntervals(data,'dist2fly_mm',...
  'nframespre',t0_on,...
  'nframespost',t1_on,...
  'on',true,...
  'nbins',nbins,...
  'binlimits',binlimits_speed);

data = HistogramFeatureIntervals(data,'dist2fly_mm',...
  'nframespre',t0_off,...
  'nframespost',t1_off,...
  'on',false,...
  'nbins',nbins,...
  'binlimits',binlimits_speed);

hfig = figure(10);
clf;
hax = PlotHistogramOnVsOff(data,'dist2fly_mm',...
  'featlabel','Inter-fly dist. (mm)');
