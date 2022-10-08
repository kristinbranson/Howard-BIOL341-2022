%% Script to post-process the output from FlyTracker and to visualize it. 

%% set up paths for code

% add JAABA to the path (.. is the parent directory to the current one)
current_directory = pwd;
parent_directory = fileparts(currentdir);
addpath(fullfile(parentDirectory, 'JAABA', 'misc'));
addpath(fullfile(parentDirectory, 'JAABA', 'filehandling'));

% add path to BIOL 341 code if not in this directory already
% uncomment this!
%addpath BIOL341;

%% set up paths for the expriment you want to look at

% where the data directories are (change)
% rootdatadir = '/groups/branson/bransonlab/alice/temp_howard';
rootdatadir = 'E:\BIOL341\GoogleDrive';
% Names of single experiments by directory name:
expnames = {
  'CsChr_JRC_SS56987_RigA_20210902T070106',
  };

% combine them to make a full path
expdirs = cell(size(expnames));
for i = 1:numel(expnames)
  expdirs{i} = fullfile(rootdatadir,expnames{i});
end
disp(expdirs)
%% look at a video

% pick a movie (1 - number of expnames)
expi = 1;
expdir = expdirs{expi};
disp(sprintf('working on exp %s', expdir))

% open playfmf with a specific video
playfmf('moviefile',fullfile(expdir,'movie.ufmf'));

% if you call just playfmf it will prompt you for the video
%playfmf;

%% post-process the tracking for further analysis

for expi = 1:numel(expdirs)
  expdir = expdirs{expi};
  fprintf('%d: %s\n',expi,expdir);
  res = PostProcessFlyTracker4JAABA(expdir);
end

%% plot an example frame with tracking

expi = 1;
fr = 100;
expdir = expdirs{expi};
data = LoadTracking(expdir);

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(fullfile(expdir,'movie.ufmf'));
im = readframe(fr);
fclose(fid);

hfig = figure(1);
clf;
hax = gca;
imagesc(im,[0,255]);
colormap gray;
hold on;
ntrajs = size(data.summary.flies,1);
colors = jet(ntrajs);
hbody = gobjects(ntrajs,1);
hwings = gobjects(ntrajs,1);
for flyi = 1:ntrajs,
  trk = data.exp.fly(flyi);
  [hbody(flyi),hwings(flyi)] = DrawWingedFly(trk.x_px(fr),trk.y_px(fr),trk.theta_rad(fr),...
    trk.a_px(fr),trk.b_px(fr),trk.xwingl_px(fr),trk.ywingl_px(fr),...
    trk.xwingr_px(fr),trk.ywingr_px(fr),'style','ellipse',...
    'Color',colors(flyi,:),'LineWidth',2);
end
axis image off;


%% plot trajectories

data = LoadTracking(expdirs);
nflies = size(data.summary.flies,1);
T = max(data.summary.exps.nframes);
nc = 5;
nr = ceil(nflies/nc);
mintimestamp = inf;
maxtimestamp = 0;
nexps = numel(data.exp);
for expi = 1:nexps,
  mintimestamp = min(mintimestamp,min(data.exp(expi).timestamps));
  maxtimestamp = max(maxtimestamp,max(data.exp(expi).timestamps));
end

hfig = figure(2);
clf;
maxr = 0;
hax = gobjects(nflies,1);
for flyi = 1:nflies,
  expnum = data.summary.flies.expnum(flyi);
  flynum = data.summary.flies.flynum(flyi);
  trk = data.exp(expnum).fly(flynum);
  hax(flyi) = subplot(nr,nc,flyi);
  plot(trk.x_mm,trk.y_mm,'k-');
  hold on;
  scatter(trk.x_mm,trk.y_mm,[],data.exp(expnum).timestamps,'.');
  title(sprintf('Exp %d, fly %d, %s',expnum,flynum,trk.sex));
  maxr = max([maxr,max(abs(trk.x_mm)),max(abs(trk.y_mm))]);
end

axis(hax,'equal');
set(hax,'XLim',[-maxr,maxr]*1.01,'YLim',[-maxr,maxr]*1.01);
set(hax,'CLim',[mintimestamp,maxtimestamp]);

%% plot speed

minspeed = 0;
maxspeed = 50; 

data = LoadTracking(expdirs);

% compute speed for each fly, and add it to the data struct.
nflies = size(data.summary.flies,1);
speed = cell(nflies,1);
for expi = 1:numel(data.exp),
  for flyi = 1:numel(data.exp(expi).fly),
    x = data.exp(expi).fly(flyi).x_mm;
    y = data.exp(expi).fly(flyi).y_mm;
    data.exp(expi).fly(flyi).speed_mmps = sqrt( (x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2 ) * data.exp(expi).summary.fps;
  end
end

hfig = figure(3);
clf;

PlotFlyFeatureOverVideo(data,'speed_mmps','minfeatplot',minspeed,'maxfeatplot',maxspeed,'featlabel','Speed (mm/s)');