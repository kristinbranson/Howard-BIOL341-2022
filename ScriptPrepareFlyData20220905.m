%% set up path

addpath /groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis;
modpath;
addpath /groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/JAABA/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/JAABA/filehandling;

%%

datafile = '/groups/branson/home/bransonk/behavioranalysis/code/MABe2022/SocialFlyBubbleExperiments_v2.csv';
cxrootdir = '/nearline/branson/FlyBubble_deep_storage';
maberootdir = '/groups/branson/home/bransonk/behavioranalysis/code/MABe2022/data';
finaldatadir = '/groups/branson/home/bransonk/behavioranalysis/code/MABe2022/sharedata20211230';
rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/BIOL341/data';


%% read in info about experiments to analyze

fid = fopen(datafile,'r');
s = fgetl(fid);
headers = regexp(s,',','split');
expinfo = [];
expdiri = find(strcmpi(headers,'file_system_path'));
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  ss = regexp(s,',','split');
  infocurr = struct;
  for j = 1:numel(ss),
    if j == expdiri,
      expdirstr = regexp(strtrim(ss{j}),'\s+','split');
      if numel(expdirstr) > 1,
        fprintf('Found spaces in file system path %s, using %s.\n',ss{j},expdirstr{1});
      end
      ss{j} = expdirstr{1};
    end
    infocurr.(headers{j}) = ss{j};
  end
  if isempty(infocurr.label),
    continue;
  end
  expinfo = structappend(expinfo,infocurr);
end

fclose(fid);

nexps = numel(expinfo);
for i = 1:nexps,
  if ~exist(expinfo(i).file_system_path,'dir'),
    fprintf('Directory %d %s missing\n',i,expinfo(i).file_system_path);
  end
end

for i = 1:nexps,
  [~,expname] = fileparts(expinfo(i).file_system_path);
  expdir = fullfile(maberootdir,expname);
  finalexpdir = fullfile(finaldatadir,expname);
  expinfo(i).isfinal = exist(finalexpdir,'dir')>0;
  expinfo(i).expdir = expdir;
  m = regexp(expname,'Rig(?<rig>.)_(?<timestamp>.*)$','names');
  expinfo(i).rig = m.rig;
  expinfo(i).timestamp = m.timestamp;
  expinfo(i).date = expinfo(i).timestamp(1:8);
end

%% which data to select

% select 2 videos per fly type from the same set
flytypes = {'aIPgBlind_newstim','pC1dpublished1_newstim','Control_RGB','BlindControl'};
idxselect = false(1,numel(expinfo));
for typei = 1:numel(flytypes),
  idx = find(strcmp({expinfo.label},flytypes{typei}) & [expinfo.isfinal]);
  fprintf('\n%s: %d\n',flytypes{typei},numel(idx));
  fprintf('  %s\n',expinfo(idx).expdir);
  ts = datenum({expinfo(idx).timestamp},'yyyymmddTHHMMSS');  
  for i = 1:numel(idx),
    idxset = find(abs(ts(i)-ts) <= 5/24/60);
    if numel(idxset) >= 2,
      idxselect(idx(idxset(1:2))) = true;
      break;
    end
  end
  fprintf('Selected:\n');
  fprintf('  %s\n',expinfo(idx(idxset(1:2))).expdir);
end

mabelabels = {expinfo(idxselect).label};
mabeexpdirs = {expinfo(idxselect).expdir};


%% back-up line

minnflies = 9;
maxnflies = 11;
maxtimediff_days = 10/24/60;

if true,
  expnames = dir(fullfile(rootdatadir,'cx_GMR_OL0046*'));
  expnames = {expnames.name};
  backupexpdirs = cellfun(@(x) fullfile(cxrootdir,x),expnames,'Uni',0);
else

  newexpdirs = mydir(fullfile(cxrootdir,'cx_GMR_OL0046*'),'isdir',true); 
  ts = nan(size(newexpdirs));
  for i = 1:numel(newexpdirs),
    [~,expname] = fileparts(newexpdirs{i});
    m = regexp(expname,'Rig._(?<timestamp>.*)$','names');
    ts(i) = datenum(m.timestamp,'yyyymmddTHHMMSS');
  end
  [ts,order] = sort(ts);
  newexpdirs = newexpdirs(order);


  ntraj = nan(size(newexpdirs));
  nflies = nan(size(newexpdirs));
  issuccess = false(size(newexpdirs));
  canchoose = false(size(newexpdirs));
  for i = 1:numel(newexpdirs),
    trxfile = fullfile(newexpdirs{i},'registered_trx.mat');
    if ~exist(trxfile,'file'),
      continue;
    end
    fprintf('%d/%d: %s\n',i,numel(newexpdirs),newexpdirs{i});
    td = load(trxfile);
    ntraj(i) = numel(td.trx);
    isalive = false([numel(td.trx),max([td.trx.endframe])]);
    for j = 1:numel(td.trx),
      isalive(j,td.trx(j).firstframe:td.trx(j).endframe) = true;
    end
    nalive = sum(isalive,1);
    nflies(i) = median(nalive);
    tmp = load(fullfile(newexpdirs{i},'automatic_checks_incoming_info.mat'));
    issuccess(i) = tmp.success>0 && exist(fullfile(newexpdirs{i},'SUCCESS'),'file');

    canchoose(i) = issuccess(i) && nflies(i) >= minnflies && nflies(i) <= maxnflies;

    idxset = find(canchoose(1:i) & abs(ts(i)-ts(1:i)) <= maxtimediff_days);
    if numel(idxset) > 1,
      newidxselect = idxset(1:2);
      break;
    end

  end

  backupexpdirs = newexpdirs(newidxselect);
end

backuplabels = repmat({'backup'},size(backupexpdirs));

allexpdirs = [mabeexpdirs,backupexpdirs];
alllabels = [mabelabels,backuplabels];

backup_analysis_protocol = '20150428_flybubble_centralcomplex';
mabe_analysis_protocol = 'current_non_olympiad_dickson_VNC';

%% copy data collection files only

copyfiles = {'movie.ufmf','Metadata.xml','Log.txt','StimulusTimingLog.txt',...
  'cx__BARCODE.csv','QuickStats.png','QuickStats.txt','SUCCESS','ABORTED',...
  'CsChr__BARCODE.csv','NorpA__BARCODE.csv','protocol.mat'};

if ~exist(rootdatadir,'dir'),
  mkdir(rootdatadir);
end

nexps = numel(allexpdirs);
expdirs = cell(1,nexps);

for i = 1:nexps,
  [~,expname] = fileparts(allexpdirs{i});
  outexpdir = fullfile(rootdatadir,expname);
  if exist(outexpdir,'dir'),
    %rmdir(outexpdir,'s');
    expdirs{i} = outexpdir;
    continue;
  end
  fprintf('Copying %d %s\n',i,allexpdirs{i});
  expdirs{i} = SymbolicCopyExperimentDirectory(allexpdirs{i},rootdatadir,'copyfiles',copyfiles);
end

%% track using flytracker

circlediameter_mm = 53.3780;
fprintf('Circle diameter should be %f\n',circlediameter_mm);

% run tracker to run tracker
tracker;

%% copy over tracking data for some of the videos
flytrackerfiles_in = {'movie-bg.mat','movie-feat.mat','movie_JAABA','movie-params.mat','movie-track.mat','flytracker-calibration.mat'};
flytrackerfiles_out = {'movie/movie-bg.mat','movie/movie-feat.mat','movie/movie_JAABA','movie/movie-params.mat','movie/movie-track.mat','calibration.mat'};
for i = 1:numel(allexpdirs),
  inexpdir = allexpdirs{i};
  outexpdir = expdirs{i};
  moviedir = fullfile(outexpdir,'movie');
  if ~exist(moviedir,'dir'),
    mkdir(moviedir);
  end
  for j = 1:numel(flytrackerfiles_in),
    infile = fullfile(inexpdir,flytrackerfiles_in{j});
    outfile = fullfile(outexpdir,flytrackerfiles_out{j});
    if exist(outfile,'file'),
      fprintf('%d,%d: %s exists, skipping\n',i,j,outfile);
      continue;
    end
    if ~exist(infile,'file'),
      fprintf('%d,%d: %s missing, skipping\n',i,j,infile);
      continue;
    end
    [success,msg] = copyfile(infile,outfile);
    if ~success,
      warning('%d,%d: failed to copy %s to %s: %s\n',i,j,infile,outfile,msg);
    end
  end
end

%% fly sizes

alla = cell(1,numel(expdirs));
allb = cell(1,numel(expdirs));
for i = 1:numel(expdirs),
  expdir = expdirs{i};
  load(fullfile(expdir,'processed_trx.mat'),'trx');
  acurr = nan(1,numel(trx));
  bcurr = nan(1,numel(trx));
  for j = 1:numel(trx),
    acurr(j) = nanmedian(trx(j).a_mm);
    bcurr(j) = nanmedian(trx(j).b_mm);
  end
  alla{i} = acurr;
  allb{i} = bcurr;
end
clf;
subplot(3,1,1);
hold on;
for i = 1:numel(expdirs),
  plot(alla{i},'.-');
end
subplot(3,1,2);
hold on;
for i = 1:numel(expdirs),
  plot(allb{i},'.-');
end
subplot(3,1,3);
hold on;
for i = 1:numel(expdirs),
  plot(alla{i}.*allb{i},'.-');
end

%% post-process data from all experiments

for i = 1:numel(expdirs),
  expdir = expdirs{i};
  postProcessFlyTracker4JAABA(expdir);
end

%% load in all data

maxT = inf;
fps = 150;

% nframesperexp = nan(size(expdirs));
% for i = 1:numel(expdirs),
%   expdir = expdirs{i};
%   headerinfo = ufmf_read_header(fullfile(expdir,'movie.ufmf'));
%   nframesperexp(i) = headerinfo.nframes;
% end

d = 5;
data = cell(1,numel(expdirs));
sex = cell(1,numel(expdirs));
names = {'x','y','ori','wing_anglel','wing_angler'};
activation = struct;
activation.startframes = cell(1,numel(expdirs));
activation.endframes = cell(1,numel(expdirs));
activation.intensities = cell(1,numel(expdirs));
activation.pulsewidths = cell(1,numel(expdirs));
activation.pulseperiods = cell(1,numel(expdirs));

for i = 1:numel(expdirs),
  expdir = expdirs{i};
  pptd = load(fullfile(expdir,'processed_trx.mat'));
  ntraj = numel(pptd.trx);

  T0 = max([pptd.trx.endframe]);
  nflies = zeros(T0,1);
  for j = 1:numel(pptd.trx),
    t0 = pptd.trx(j).firstframe;
    t1 = pptd.trx(j).endframe;
    nflies(t0:t1) = nflies(t0:t1)+1;
  end
  mednflies = median(nflies);
  T = find(nflies>=mednflies,1,'last');
  T = min(maxT,T);
  fprintf('%d: T0 = %d, T = %d\n',i,T0,T);

  datacurr = nan([T,d,ntraj]);
  sexcurr = repmat('?',[1,ntraj]);
  load(fullfile(expdir,'indicatordata.mat'),'indicatorLED');
  startidx = indicatorLED.startframe <= T;
  endidx = indicatorLED.endframe <= T;
  assert(nnz(startidx)==nnz(endidx));
  idx = startidx;
  activation.startframes{i} = indicatorLED.startframe(idx);
  activation.endframes{i} = indicatorLED.endframe(idx);
  activation.intensities{i} = indicatorLED.intensity(idx);
  activation.pulsewidths{i} = indicatorLED.pulsewidths(idx);
  activation.pulseperiods{i} = indicatorLED.pulseperiods(idx);

  for j = 1:numel(pptd.trx),
    x = nan(T,1);
    y = nan(T,1);
    theta = nan(T,1);
    wing_anglel = nan(T,1);
    wing_angler = nan(T,1);
    t0 = pptd.trx(j).firstframe;
    t1 = min(pptd.trx(j).endframe,T);
    i1 = t1-t0+1;
    x(t0:t1) = pptd.trx(j).x_mm(1:i1);
    y(t0:t1) = pptd.trx(j).y_mm(1:i1);
    theta(t0:t1) = pptd.trx(j).theta_mm(1:i1);
    wing_anglel(t0:t1) = pptd.trx(j).wing_anglel(1:i1);
    wing_angler(t0:t1) = pptd.trx(j).wing_anglel(1:i1);
    datacurr(:,:,j) = [x,y,theta,wing_anglel,wing_angler];
    sexcurr(j) = pptd.trx(j).sex;
  end
  data{i} = datacurr;
  sex{i} = sexcurr;
end
expnames = cell(1,numel(expdirs));
for i = 1:numel(expnames),
  [~,expnames{i}] = fileparts(expdirs{i});
end
%save('SampleData20220912.mat','data','names','sex','expnames','activation','fps');

%% select intervals 

pulsefrac = 1;
actintensity = 'max';
startofftime = 10*fps;

datastartframes = nan(1,numel(expdirs));
dataendframes = nan(1,numel(expdirs));
for i = 1:numel(expdirs),

  non = numel(activation.intensities{i});
  idxon = true(1,non);
  if ~isempty(pulsefrac),
    idxon = (activation.pulsewidths{i} ./ activation.pulseperiods{i}) == pulsefrac;
  end
  switch actintensity,
    case 'max',
      [maxint] = max(activation.intensities{i}(idxon));
      idxon = idxon & activation.intensities{i} == maxint;
    case 'all'
      % do nothing
  end

  idxon = find(idxon);
  if idxon(1) == 1,
    t0 = 1;
  else
    t0 = activation.endframes{i}(idxon(1)-1)+1+startofftime;
  end
  if idxon(end) == non,
    t1 = size(data{i},1);
  else 
    t1 = activation.startframes{i}(idxon(end)+1)-1;
  end

  datastartframes(i) = t0;
  dataendframes(i) = t1;
  data{i} = data{i}(t0:t1,:,:);
  activation.intensities{i} = activation.intensities{i}(idxon);
  activation.pulsewidths{i} = activation.pulsewidths{i}(idxon);
  activation.pulseperiods{i} = activation.pulseperiods{i}(idxon);
  activation.startframes{i} = activation.startframes{i}(idxon)-t0+1;
  activation.endframes{i} = activation.endframes{i}(idxon)-t0+1;

end

%% save to file

save('SampleData20220915.mat','data','names','sex','expnames','activation','fps','datastartframes','dataendframes');