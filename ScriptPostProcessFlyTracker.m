%% Script to post-process the output from FlyTracker and to visualize it. 

%% set up paths for code

% add JAABA to the path (assume it is in the parent directory to the current one)
current_directory = pwd;
parent_directory = fileparts(current_directory);
addpath(fullfile(parent_directory, 'JAABA', 'misc'));
addpath(fullfile(parent_directory, 'JAABA', 'filehandling'));

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
fprintf('Here is the first full exp dir:\n%s\n', expdirs{i})
%% look at a video

% pick a movie (1 - number of expnames)
expi = 1;
expdir = expdirs{expi};
fprintf('working on exp %s\n', expdir)

% open playfmf with a specific video (assume it is called movie.ufmf)
movie_name = fullfile(expdir,'movie.ufmf');
fprintf('trying to open movie:\n %s\n', movie_name)
playfmf('moviefile',movie_name);

% if you call just playfmf it will prompt you for the video
%playfmf;

%% post-process the tracking for further analysis

for expi = 1:numel(expdirs)
  expdir = expdirs{expi};
  fprintf('Working on #%d: %s\n',expi,expdir);
  res = PostProcessFlyTracker4JAABA(expdir);
end

%% plot an example frame with tracking

% set up name and frame to look at
expi = 1;
fr = 107;
expdir = expdirs{expi};
fprintf('working on exp %s\n', expdir)
movie_name = fullfile(expdir,'movie.ufmf');
fprintf('trying to open frame #%d from file movie:\n%s\n', fr, movie_name)

% get a function (readframe) according to the movie you want to read
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(movie_name);
% read the frame
im = readframe(fr);
sz = size(im);
fprintf('Read 1 frame size is %d x %d pixels \n', sz(1), sz(2))
% close the file (no more reading)
fclose(fid);

% plot the frame in grayscale
figure(1);
clf;
imagesc(im,[0,255]);
colormap gray;

% add tracking information, what is LoadTracking? mark the function and
% press F1 to get help. Also look at the result 'data' in the Variable
% explorer

data = LoadTracking(expdir); % explor this struct in the workspace!
hold on;
ntrajs = size(data.summary.flies,1);
fprintf('Found %d files in the metadata\n', ntrajs);

colors = jet(ntrajs); % get a color for each one
% these next 2 lines make array in a type that can hold graphic objects
% (not arrays of numbers like zeros or cell arrays)
hbody = gobjects(ntrajs,1); 
hwings = gobjects(ntrajs,1);
% loop over each track coresponding to a fly
for flyi = 1:ntrajs
    % get the track
    trk = data.exp.fly(flyi);
    % get locations based on the frame we want (fr)
    x = trk.x_px(fr);
    y = trk.y_px(fr);
    fprintf('Drawing fly #%d at location x: %1.1f, y: %1.1f in pixels\n', flyi, x, y)
    angle = trk.theta_rad(fr);
    a =  trk.a_px(fr);
    b = trk.b_px(fr);
    x_wingl = trk.xwingl_px(fr);
    y_wingl = trk.ywingl_px(fr);
    xwingr = trk.xwingr_px(fr);
    ywingr = trk.ywingr_px(fr);
    color = colors(flyi,:);
    % use all that data for the function DrawWingedFly that returns the
    [hbody(flyi),hwings(flyi)] = DrawWingedFly( ...
        x, y, angle, a, b, ...
        x_wingl, y_wingl, xwingr, ywingr, ...
        'style','ellipse', 'Color',color,'LineWidth',2);
end
axis image off;
% worked until here, I have only one example movie so 

%% plot trajectories

data = LoadTracking(expdirs);
nflies = size(data.summary.flies,1);
T = max(data.summary.exps.nframes); % max number of frames in the movie
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