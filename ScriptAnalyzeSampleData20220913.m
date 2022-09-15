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
t0 = 1; % start frame to plot
nframesplot = 3*60*fps; % number of frames to plot
videoi = 1; % which video to plot
fly = 8; % which fly to plot

hfig = figure(1);
clf;
t1 = min(size(data{videoi},1),t0 + nframesplot - 1);
plot(data{videoi}(t0:t1,XIDX,fly),data{videoi}(t0:t1,YIDX,fly),'k-');
hold on;
scatter(data{videoi}(t0:t1,XIDX,fly),data{videoi}(t0:t1,YIDX,fly),[],(t0:t1)/fps,'.');
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
  % create the axes
  hax(i) = subplot(nfeat,1,i);
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

% average over flies
meanspeedpervideo = cell(1,nvideos);
for i = 1:nvideos,
  meanspeedpervideo{i} = mean(speed{i},2,'omitnan');
end

% compute the standard deviation over all flies
stdspeedpervideo = cell(1,nvideos);
for i = 1:nvideos,
  stdspeedpervideo{i} = std(speed{i},0,2,'omitnan');
end

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

% compute the average over all flies
meandist2flypervideo = cell(1,nvideos);
for i = 1:nvideos,
  meandist2flypervideo{i} = mean(dist2fly{i},2,'omitnan');
end

% compute the standard deviation over all flies
stddist2flypervideo = cell(1,nvideos);
for i = 1:nvideos,
  stddist2flypervideo{i} = std(dist2fly{i},0,2,'omitnan');
end


%% plot the flies' speeds for each video

hfig = figure(3);
clf;
naxc = 2;
naxr = ceil(nvideos/naxc);
maxspeedplot = 70; % limit of the y-axis in mm/s
lightcolor = [1,.7,.7]; % pink
plotallflies = false; % whether to plot individual flies
plotstderr = ~plotallflies; % whether to plot the standard deviation

% one set of axes per video
hax = gobjects(nvideos,1);
for i = 1:nvideos,
  % create the axes
  hax(i) = subplot(naxr,naxc,i);
  hold(hax(i),'on');

  % plot when the lights are on
  for j = 1:numel(activation.startframes{i}),
    t0 = activation.startframes{i}(j);
    t1 = activation.endframes{i}(j);
    patch(hax(i),[t0,t0,t1,t1,t0]/fps,[0,1,1,0,0]*maxspeedplot,lightcolor,'LineStyle','none');
  end

  % plot the speed for all the flies. this will create a different line for
  % each fly
  T = size(data{i},1);
  if plotallflies,
    plot(hax(i),(1:T-1)/fps,speed{i}','-');
    statcolor = 'k';
  else
    statcolor = genotypecolors(genotypeidx(i),:);
  end

  % plot standard error (standard deviation of mean)
  if plotstderr,
    stderrcurr = stdspeedpervideo{i}/sqrt(ntrajs(i));
    plot(hax(i),(1:T-1)/fps,meanspeedpervideo{i}-stderrcurr,'-','Color',statcolor);
    plot(hax(i),(1:T-1)/fps,meanspeedpervideo{i}+stderrcurr,'-','Color',statcolor);
  end

  % plot starts of activation periods
  for j = 1:numel(activation.startframes{i}),
    t0 = activation.startframes{i}(j);
    plot(hax(i),[t0,t0]/fps,[0,maxspeedplot],'r-','LineWidth',2);
  end

  % plot the mean speed over all flies
  plot(hax(i),(1:T-1)/fps,meanspeedpervideo{i}','-','Linewidth',2,'Color',statcolor);

  title(hax(i),sprintf('Video %d (%s)',i,expnames{i}),'Interpreter','none');
  box(hax(i),'off');
end
% all axes forced to have the same x and y limits
linkaxes(hax);
xlabel(hax(end),'Time (s)');
ylabel(hax(end),'Speed (mm/s)');
set(hax,'YLim',[0,maxspeedplot],'XLim',[0,(T+1)/fps])

%% plot the distance apart for each fly over the videos

hfig = figure(4);
clf;
mindistplot = 0; % limits of the y-axis in mm
maxdistplot = 15;
plotallflies = false; % whether to plot individual flies
plotstderr = ~plotallflies; % whether to plot the standard deviation
lightcolor = [1,.7,.7]; % pink

% one set of axes per video
hax = gobjects(nvideos,1);
for i = 1:nvideos,
  % create the axes
  hax(i) = subplot(nvideos,1,i);
  hold(hax(i),'on');

  % plot when the lights are on
  for j = 1:numel(activation.startframes{i}),
    t0 = activation.startframes{i}(j);
    t1 = activation.endframes{i}(j);
    patch(hax(i),[t0,t0,t1,t1,t0]/fps,...
      [mindistplot,maxdistplot,maxdistplot,mindistplot,mindistplot],...
      lightcolor,'LineStyle','none');
  end

  % plot the inter-fly distance for all the flies. this will create a
  % different line for each fly
  if plotallflies,
    plot(hax(i),(1:T)/fps,dist2fly{i}','-');
    statcolor = 'k';
  else
    statcolor = genotypecolors(genotypeidx(i),:);
  end
  % plot standard error (standard deviation of mean)
  if plotstderr,
    stderrcurr = stddist2flypervideo{i}/sqrt(ntrajs(i));
    plot(hax(i),(1:T)/fps,meandist2flypervideo{i}-stderrcurr,'-','Color',statcolor);
    plot(hax(i),(1:T)/fps,meandist2flypervideo{i}+stderrcurr,'-','Color',statcolor);
  end

  % plot starts of activation periods
  for j = 1:numel(activation.startframes{i}),
    t0 = activation.startframes{i}(j);
    plot(hax(i),[t0,t0]/fps,[mindistplot,maxdistplot],'r-','LineWidth',2);
  end

  % plot the mean inter-fly distance over all flies
  plot(hax(i),(1:T)/fps,meandist2flypervideo{i}','-','Linewidth',2,'Color',statcolor);

  title(hax(i),sprintf('Video %d (%s)',i,expnames{i}),'Interpreter','none');
  box(hax(i),'off');
end

% all axes forced to have the same x and y limits
linkaxes(hax);
xlabel(hax(end),'Time (s)');
ylabel(hax(end),'Inter-fly distance (mm)');
set(hax,'YLim',[mindistplot,maxdistplot],'XLim',[0,(T+1)/fps])

%% plot the speeds when the LED turns on

nframespre_plot = 1*fps; % how much before the lights on to plot
nframespost_plot = 5*fps; % how much after the lights on to plot
nstimuliplot = 13; % how many stimulus periods to plot, should be a number between 1 and 13
maxspeedplot = 40; % maximum speed to plot
plotallflies = false; % whether to plot individual flies
plotstderr = ~plotallflies; % whether to plot the standard deviation

hfig = figure(5);
clf;
% row of axes for each video, column for each stimulus period
hax = gobjects(nvideos,nstimuliplot);

for i = 1:nvideos,
  for j = 1:nstimuliplot,
    if j > numel(activation.startframes{i}),
      continue;
    end

    % which axes to plot in
    axi = sub2ind([nstimuliplot,nvideos],j,i);
    hax(i,j) = subplot(nvideos,nstimuliplot,axi);

    % which frames to plot
    t = activation.startframes{i}(j);
    t0 = max(1,t-nframespre_plot);
    t1 = min(T,t+nframespost_plot);

    hold(hax(i,j),'on');
    % plot each fly
    if plotallflies,
      plot(hax(i,j),(t0-t:t1-t)/fps,speed{i}(t0:t1,:)','-');
      statcolor = 'k';
    else
      statcolor = genotypecolors(genotypeidx(i),:);
    end
    if plotstderr,
      stderrcurr = stdspeedpervideo{i}(t0:t1)/sqrt(ntrajs(i));
      plot(hax(i,j),(t0-t:t1-t)/fps,meanspeedpervideo{i}(t0:t1)-stderrcurr,'-','Color',statcolor);
      plot(hax(i,j),(t0-t:t1-t)/fps,meanspeedpervideo{i}(t0:t1)+stderrcurr,'-','Color',statcolor);
    end

    % plot a vertical line when the lights turn on
    plot(hax(i,j),[0,0],[0,maxspeedplot],'r-');
    % plot the average speed of all flies
    if plotallflies,
      lw = 2;
    else
      lw = 1;
    end
    plot(hax(i,j),(t0-t:t1-t)/fps,meanspeedpervideo{i}(t0:t1)','-','LineWidth',lw,'color',statcolor);

    % only show which stimulus this is in the top row of axes
    if i == 1,
      title(hax(i,j),sprintf('Stimulus %d',j));
    end
    box(hax(i,j),'off');
  end
  % only show which video this is in the left column of axes
  ylabel(hax(i,1),sprintf('Video %d',i));
end

% force axes to share the same x and y limits
linkaxes(hax(ishandle(hax)));
ylabel(hax(end,1),{'Speed (mm/s)',sprintf('Video %d',nvideos)});
xlabel(hax(end,1),'Time (s)');
set(hax(ishandle(hax)),'YLim',[0,maxspeedplot],'XLim',[-nframespre_plot,nframespost_plot]/fps);


%% plot the inter-fly distances when the LED turns on

nframespre_plot = 5*fps; % how much before the lights on to plot
nframespost_plot = 30*fps; % how much after the lights on to plot
nstimuliplot = 13; % how many stimulus periods to plot, should be a number between 1 and 13
mindistplot = 0; % limits of the y-axis in mm
maxdistplot = 15;
plotallflies = false; % whether to plot individual flies
plotstderr = ~plotallflies; % whether to plot the standard deviation

hfig = figure(6);
clf;

% row of axes for each video, column for each stimulus period
hax = gobjects(nvideos,nstimuliplot);

for i = 1:nvideos,
  for j = 1:nstimuliplot,
    if j > numel(activation.startframes{i}),
      continue;
    end
    % which axes to plot in
    axi = sub2ind([nstimuliplot,nvideos],j,i);
    hax(i,j) = subplot(nvideos,nstimuliplot,axi);

    % which frames to plot
    t = activation.startframes{i}(j);
    t0 = max(1,t-nframespre_plot);
    t1 = min(T,t+nframespost_plot);
    
    hold(hax(i,j),'on');

    % plot a line for each fly
    if plotallflies,
      plot(hax(i,j),(t0-t:t1-t)/fps,dist2fly{i}(t0:t1,:)','-');
      statcolor = 'k';
    else
      statcolor = genotypecolors(genotypeidx(i),:);
    end
    % plot standard error
    if plotstderr,
      stderrcurr = stddist2flypervideo{i}(t0:t1)/sqrt(ntrajs(i));
      plot(hax(i,j),(t0-t:t1-t)/fps,meandist2flypervideo{i}(t0:t1)-stderrcurr,'-','Color',statcolor);
      plot(hax(i,j),(t0-t:t1-t)/fps,meandist2flypervideo{i}(t0:t1)+stderrcurr,'-','Color',statcolor);
    end

    % plot a vertical line indicating lights on time
    plot(hax(i,j),[0,0],[mindistplot,maxdistplot],'r-');

    % plot average over flies
    if plotallflies,
      lw = 2;
    else
      lw = 1;
    end
    plot(hax(i,j),(t0-t:t1-t)/fps,meandist2flypervideo{i}(t0:t1)','-','LineWidth',lw,'Color',statcolor);

    % only show which stimulus this is in the top row of axes
    if i == 1,
      title(hax(i,j),sprintf('Stimulus %d',j));
    end
    box(hax(i,j),'off');
  end
  % only show which video this is in the left column of axes
  ylabel(hax(i,1),sprintf('Video %d',i));
end

% force all axes to share the same x and y limits
linkaxes(hax(ishandle(hax)));

ylabel(hax(end,1),{'Inter-fly distance (mm)',sprintf('Video %d',nvideos)});
xlabel(hax(end,1),'Time (s)');
set(hax(ishandle(hax)),'YLim',[mindistplot,maxdistplot],'XLim',[-nframespre_plot,nframespost_plot]/fps);

%% plot the average speed over all stimuli

nframespre_plot = 1*fps; % how much before the lights on to plot
nframespost_plot = 5*fps; % how much after the lights on to plot
maxspeedplot = 25; % maximum speed to plot
plotstderr = true; % whether to plot the standard error

hfig = figure(7);
clf;

hax = gobjects(nvideos,1);
for i = 1:nvideos,

  videocolor = genotypecolors(genotypeidx(i),:);
  meanspeedon = zeros(nframespre_plot+nframespost_plot+1,1);
  count = zeros(nframespre_plot+nframespost_plot+1,1);
  for j = 1:numel(activation.startframes{i}),
    t = activation.startframes{i}(j);
    t0 = max(1,t-nframespre_plot);
    t1 = min(T,t+nframespost_plot);
    i0 = t0-(t-nframespre_plot) + 1;
    i1 = i0 + (t1-t0); 
    meanspeedon(i0:i1) = meanspeedon(i0:i1) + sum(speed{i}(t0:t1,:),2,'omitnan');
    count(i0:i1) = count(i0:i1) + sum(~isnan(speed{i}(t0:t1,:)),2);
  end
  meanspeedon = meanspeedon ./ count;
  stdspeedon = zeros(nframespre_plot+nframespost_plot+1,1);
  for j = 1:numel(activation.startframes{i}),
    t = activation.startframes{i}(j);
    t0 = max(1,t-nframespre_plot);
    t1 = min(T,t+nframespost_plot);
    i0 = t0-(t-nframespre_plot) + 1;
    i1 = i0 + (t1-t0); 
    stdspeedon(i0:i1) = stdspeedon(i0:i1) + sum((speed{i}(t0:t1,:)-meanspeedon(i0:i1)).^2,2,'omitnan');
  end
  stdspeedon = sqrt(stdspeedon ./ count);
  hax(i) = subplot(nvideos,1,i);
  plot([0,0],[0,maxspeedplot],'k-');
  hold on;

  if plotstderr,
    plot((-nframespre_plot:nframespost_plot)/fps,meanspeedon-stdspeedon./max(1,sqrt(count)),'-','Color',videocolor);
    plot((-nframespre_plot:nframespost_plot)/fps,meanspeedon+stdspeedon./max(1,sqrt(count)),'-','Color',videocolor);
  end
  plot((-nframespre_plot:nframespost_plot)/fps,meanspeedon,'-','Color',videocolor,'LineWidth',2);
  title(sprintf('Video %d: %s',i,expnames{i}),'Interpreter','none');
end
ylabel(hax(end),'Speed (mm/s)');
xlabel(hax(end),'Time (s)');
set(hax,'YLim',[0,maxspeedplot],'XLim',[-nframespre_plot,nframespost_plot]/fps);

%% plot the average inter fly distance over all stimuli

nframespre_plot = 5*fps; % how much before the lights on to plot
nframespost_plot = 30*fps; % how much after the lights on to plot
nstimuliplot = 13; % how many stimulus periods to plot, should be a number between 1 and 13
mindistplot = 0; % limits of the y-axis in mm
maxdistplot = 15;
plotstderr = true; % whether to plot the standard error

hfig = figure(8);
clf;

hax = gobjects(nvideos,1);
for i = 1:nvideos,

  videocolor = genotypecolors(genotypeidx(i),:);
  meandist2flyon = zeros(nframespre_plot+nframespost_plot+1,1);
  count = zeros(nframespre_plot+nframespost_plot+1,1);
  for j = 1:numel(activation.startframes{i}),
    t = activation.startframes{i}(j);
    t0 = max(1,t-nframespre_plot);
    t1 = min(T,t+nframespost_plot);
    i0 = t0-(t-nframespre_plot) + 1;
    i1 = i0 + (t1-t0); 
    meandist2flyon(i0:i1) = meandist2flyon(i0:i1) + sum(dist2fly{i}(t0:t1,:),2,'omitnan');
    count(i0:i1) = count(i0:i1) + sum(~isnan(dist2fly{i}(t0:t1,:)),2);
  end
  meandist2flyon = meandist2flyon ./ count;
  stddist2flyon = zeros(nframespre_plot+nframespost_plot+1,1);
  for j = 1:numel(activation.startframes{i}),
    t = activation.startframes{i}(j);
    t0 = max(1,t-nframespre_plot);
    t1 = min(T,t+nframespost_plot);
    i0 = t0-(t-nframespre_plot) + 1;
    i1 = i0 + (t1-t0); 
    stddist2flyon(i0:i1) = stddist2flyon(i0:i1) + sum((dist2fly{i}(t0:t1,:)-meandist2flyon(i0:i1)).^2,2,'omitnan');
  end
  stddist2flyon = sqrt(stddist2flyon ./ count);
  hax(i) = subplot(nvideos,1,i);
  plot([0,0],[mindistplot,maxdistplot],'k-');
  hold on;

  if plotstderr,
    plot((-nframespre_plot:nframespost_plot)/fps,meandist2flyon-stddist2flyon./max(1,sqrt(count)),'-','Color',videocolor);
    plot((-nframespre_plot:nframespost_plot)/fps,meandist2flyon+stddist2flyon./max(1,sqrt(count)),'-','Color',videocolor);
  end
  plot((-nframespre_plot:nframespost_plot)/fps,meandist2flyon,'-','Color',videocolor,'LineWidth',2);
  title(sprintf('Video %d: %s',i,expnames{i}),'Interpreter','none');
end
ylabel(hax(end),'Inter-fly distance (mm)');
xlabel(hax(end),'Time (s)');
set(hax,'YLim',[mindistplot,maxdistplot],'XLim',[-nframespre_plot,nframespost_plot]/fps);

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

  % histogram each on interval for each fly
  for j = 1:numel(onstarts),
    % which frames to histogram
    t = onstarts(j);
    t0 = min(T-1,max(1,t+t0_on));
    t1 = min(T-1,max(1,t+t1_on));

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
    t0 = min(T-1,max(1,t+t0_off));
    t1 = min(T-1,max(1,t+t1_off));

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