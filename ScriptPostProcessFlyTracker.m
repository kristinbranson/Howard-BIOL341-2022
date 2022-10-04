%% Script to post-process the output from FlyTracker and to visualize it. 

%% set up paths

% add JAABA to the path
addpath ../JAABA/misc;
addpath ../JAABA/filehandling;

% add path to BIOL 341 code if not in this directory already
% uncomment this!

%addpath BIOL341;

rootdatadir = '/groups/branson/bransonlab/alice/temp_howard';

% set data locations
expnames = {
  'HU_Back_Ctrl_RigF_20220927T172138',
  'HU_Back_Ctrl_RigE_20220927T172346',
  'CsChrSocial3_aIPg_RigF_20220923T163615'
  };

expdirs = cell(size(expnames));
for i = 1:numel(expnames),
  expdirs{i} = fullfile(rootdatadir,expnames{i});
end

%% look at a video

expi = 1;
expdir = expdirs{expi};

% open playfmf with a specific video
playfmf('moviefile',fullfile(expdir,'movie.ufmf'));

% if you call just playfmf it will prompt you for the video
%playfmf;

%% post-process the tracking for further analysis

for expi = 1:numel(expdirs),
  expdir = expdirs{expi};
  fprintf('%d: %s\n',expi,expdir);
  res = PostProcessFlyTracker4JAABA(expdir);
end

%% plot an example frame with tracking

expi = 1;
fr = 100;
expdir = expdirs{expi};
trx = LoadTracking(expdir);

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(fullfile(expdir,'movie.ufmf'));
im = readframe(fr);
fclose(fid);

hfig = figure(1);
clf;
hax = gca;
imagesc(im,[0,255]);
colormap gray;
hold on;
colors = jet(size(trx,1));
hbody = gobjects(size(trx,1),1);
hwings = gobjects(size(trx,1),1);
for fly = 1:size(trx,1),
  [hbody(fly),hwings(fly)] = DrawWingedFly(trx.x_px(fly,fr),trx.y_px(fly,fr),trx.theta_rad(fly,fr),...
    trx.a_px(fly,fr),trx.b_px(fly,fr),trx.xwingl_px(fly,fr),trx.ywingl_px(fly,fr),...
    trx.xwingr_px(fly,fr),trx.ywingr_px(fly,fr),'style','ellipse',...
    'Color',colors(fly,:));
end
axis image off;


%% plot trajectories

trx = LoadTracking(expdirs);
[nflies,T] = size(trx.x_mm);
nc = 5;
nr = ceil(nflies/nc);
mintimestamp = min(trx.timestamp(:));
maxtimestamp = max(trx.timestamp(:));

hfig = figure(2);
clf;
maxr = max(max(abs(trx.x_mm(:))),max(abs(trx.y_mm(:))));
hax = gobjects(nflies,1);
for fly = 1:nflies,
  hax(fly) = subplot(nr,nc,fly);
  plot(trx.x_mm(fly,:),trx.y_mm(fly,:),'k-');
  hold on;
  scatter(trx.x_mm(fly,:),trx.y_mm(fly,:),[],trx.timestamp(fly,:),'.');
  set(gca,'XLim',[-maxr,maxr]*1.01,'YLim',[-maxr,maxr]*1.01);
  axis('equal');
  title(sprintf('Exp %d, fly %d, %s',trx.metadata.exp_num(fly),trx.metadata.id(fly),trx.metadata.sex(fly)));
end

set(hax,'CLim',[mintimestamp,maxtimestamp]);

%% plot speed

minspeed = 0;
maxspeed = 3; 
speed = sqrt((trx.x_mm(:,2:end)-trx.x_mm(:,1:end-1)).^2 + (trx.y_mm(:,2:end)-trx.y_mm(:,1:end-1)).^2);

hfig = figure(3);
clf;

PlotFlyFeatureOverVideo(speed,trx.metadata.fps(1),trx.metadata.activation_startframes,trx.metadata.activation_endframes,...
  'minfeatplot',minspeed,'maxfeatplot',maxspeed,'featlabel','Speed (mm/s)');