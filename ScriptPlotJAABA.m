%% set some parameters

% set where your jab file is
jabfile = 'sample_processed_data\fencing.jab';

% set parameters
t0_on = round(0); % how long after lights on to include
t1_on = inf; % maximum time after lights on to include. set to inf to include to the end of the period. 
t0_off = round(.25*30); % how long after lights off to include
t1_off = inf; % maximum time after lights off to include. set to inf to include to the end of the period. 

%% determine experiment directories and scores file
% load in data

jd = load(jabfile,'-mat');
expdirs = jd.x.expDirNames;
nexps = numel(expdirs);
scorefilename = jd.x.file.scorefilename{1};
behaviorname = jd.x.behaviors.names{1};

% make sure scores are available for all experiment directories
for expnum = 1:nexps,
  scoresfile = fullfile(expdirs{expnum},scorefilename);
  assert(exist(scoresfile,'file')>0,sprintf('scoresfile %s does not exist',scoresfile))
end

data = LoadTracking(expdirs);
data = LoadJAABAScores(data,expdirs,scorefilename);
fps = nanmean(data.summary.exps.fps);

fprintf('behavior: %s, scorefile: %s, %d experiments:\n',behaviorname,...
  scorefilename,nexps);
fprintf('%s\n',expdirs{:})

%% colors for plotting

exptypes = data.summary.exps.type;
[unique_exptypes,~,exptypeidx] = unique(exptypes);

exptypecolors = lines(numel(exptypes));

%% plot an "ethogram"

expnum = 1;
nflies = numel(data.exp(expnum).fly);
nframes = data.exp(expnum).summary.nframes;
isbehavior = zeros(nflies,nframes);
for flynum = 1:nflies,
  pred = data.exp(expnum).fly(flynum).(behaviorname);
  isbehavior(flynum,1:numel(pred)) = pred;
end

hfig = 101;
figure(hfig);
clf;
imagesc(isbehavior,'xdata',[1/fps,nframes/fps]);
colormap gray;
hold on;

% plot when the lights are on
ylim = [.5,nflies+.5];
h = PlotActivationTimePatch(data.exp(expnum).activation.startframe,data.exp(expnum).activation.endframe,fps,ylim,gca);
set(h,'FaceColor',[1,0,0],'FaceAlpha',.3);
ylabels = repmat("",nflies,1);
for flynum = 1:nflies,
  ylabels(flynum) = sprintf("%d, %s",flynum,data.exp(expnum).fly(flynum).sex);
end
set(gca,'YTick',1:nflies,'YTickLabel',ylabels);
title(sprintf('%s ethogram, %s',behaviorname,data.summary.exps.type(expnum)),'Interpreter','none');
xlabel('Time (s)');

%% compute the fraction of time performing each behavior

for expnum = 1:nexps,
  nflies = numel(data.exp(expnum).fly);

  % keep a total of the number of frames of the behavior over all flies in
  % the experiment, females and males
  countpos_any = [0,0];
  total_any = [0,0];
  countpos_on = [0,0];
  total_on = [0,0];
  countpos_off = [0,0];
  total_off = [0,0];


  % loop over flies
  for flynum = 1:nflies,

    % predictions for this behavior
    pred = data.exp(expnum).fly(flynum).(behaviorname);
    sex = double(data.exp(expnum).fly(flynum).sex == 'm')+1;

    % compute the fraction of time over all frames
    [data.exp(expnum).fly(flynum).(['fractime_',behaviorname]),...
      countpos_any_curr,total_any_curr] = ...
    ComputeFractimeActivation(pred,1,numel(pred),0,0);
    countpos_any(sex) = countpos_any(sex) + countpos_any_curr;
    total_any(sex) = total_any(sex) + total_any_curr;

    % compute the fraction of time during activation
    t0s = data.exp(expnum).activation.startframe;
    t1s = data.exp(expnum).activation.endframe;
    [data.exp(expnum).fly(flynum).(['fractime_on_',behaviorname]),...
      countpos_on_curr,total_on_curr] = ...
      ComputeFractimeActivation(pred,t0s,t1s,t0_on,t1_on);
    countpos_on(sex) = countpos_on(sex) + countpos_on_curr;
    total_on(sex) = total_on(sex) + total_on_curr;

    % compute the fraction of time not during activation
    t0s = [1;data.exp(expnum).activation.endframe+1];
    t1s = [data.exp(expnum).activation.startframe;numel(pred)];
    [data.exp(expnum).fly(flynum).(['fractime_off_',behaviorname]),...
      countpos_off_curr,total_off_curr] = ...
      ComputeFractimeActivation(pred,t0s,t1s,t0_off,t1_off);
    countpos_off(sex) = countpos_off(sex) + countpos_off_curr;
    total_off(sex) = total_off(sex) + total_off_curr;

  end

  % compute fraction time for all flies in the experiment
  data.exp(expnum).stat.(['fractime_',behaviorname]) = sum(countpos_any) / sum(total_any);
  data.exp(expnum).stat.(['fractime_on_',behaviorname]) = sum(countpos_on) / sum(total_on);
  data.exp(expnum).stat.(['fractime_off_',behaviorname]) = sum(countpos_off) / sum(total_off);

  % female
  sex = 1;
  data.exp(expnum).stat.(['fractime_',behaviorname,'_female']) = countpos_any(sex) / total_any(sex);
  data.exp(expnum).stat.(['fractime_on_',behaviorname,'_female']) = countpos_on(sex) / total_on(sex);
  data.exp(expnum).stat.(['fractime_off_',behaviorname,'_female']) = countpos_off(sex) / total_off(sex);

  % male
  sex = 2;
  data.exp(expnum).stat.(['fractime_',behaviorname,'_male']) = countpos_any(sex) / total_any(sex);
  data.exp(expnum).stat.(['fractime_on_',behaviorname,'_male']) = countpos_on(sex) / total_on(sex);
  data.exp(expnum).stat.(['fractime_off_',behaviorname,'_male']) = countpos_off(sex) / total_off(sex);

end

%% plot the fraction of time the behavior occurs

hfig = 102;
figure(hfig);
clf;

if ~isempty(data.summary.activation),
  hax = subplot(3,1,1);
else
  hax = gca;
end

% plot fraction of time for all frames
fn = ['fractime_',behaviorname];
[~,noise] = PlotFracTime(data,fn,'ax',hax);
title(hax,sprintf('%s, all frames',behaviorname));

% plot fraction of time for on frames
if ~isempty(data.summary.activation),
  hax(2) = subplot(3,1,2);
  fn = ['fractime_on_',behaviorname];
  PlotFracTime(data,fn,'ax',hax(2),'noise',noise);
  title(hax,sprintf('%s, activation on',behaviorname));

  % plot fraction of time for off frames
  hax(3) = subplot(3,1,3);
  fn = ['fractime_off_',behaviorname];
  PlotFracTime(data,fn,'ax',hax(3),'noise',noise);
  title(hax,sprintf('%s, activation off',behaviorname));
end

linkaxes(hax);
