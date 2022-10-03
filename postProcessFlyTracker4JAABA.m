% res = postProcessFlyTracker4JAABA(expdir,...)
%
% Processes trajectories and metadata in the input experiment directory:
% 1. Removes nans that can affect downstream processing. For short
% sequences of nans (<= maxFlyTrackerNanInterpFrames), it interpolates. For
% longer sequences of nans, it breaks trajectories into two tracklets. 
% 2. Classifies fly sex. If the metadata says that this is a single-sex
% video, it uses that sex. If it is a mixed sex video, uses
% FlyDiscoClassifySex to classify sex based on area of flies. 
% Processed trajectories are save to <expdir>/processed_trx.mat and
% <expdir>/processed-movie-track.mat. 
% 3. Computes wing features used by JAABA. Outputs are saved in the
% directory <expdir>/perframe. 
% 4. Determines onset and offset of activation. Outputs are saved in
% <expdir>/indicatordata.mat. 
%
% Inputs:
% expdir: path to experiment directory
%
% Output:
% res: struct with some of the outputs. All of the outputs are saved to
% file as described aboved.
%
% Optional inputs:
% nmale: For mixed sex experiments, if you know how many male flies are
% there are, you can specify this. Then, this the smallest nmale flies will
% be assigned sex 'm' and the largest 'f'. Otherwise, the clustering
% algorithm in FlyDiscoClassifySex will be used. Default: [], meaning
% unspecified. 
% Advanced optional inputs (you should not need to change these). 
% maxFlyTrackerNanInterpFrames: Maximum number of nan frames to interpolate
% through. Default: 5. 
% The following define the names of files within the experiment directory:
% trxfilestr: where FlyTracker stored the JAABA-compatible trajectories.
% Default: 'movie/movie_JAABA/trx.mat'.
% flytrackerstr: where FlyTracker stored the trajectories.
% Default: 'movie/movie-track.mat'. 
% calibfilestr: where FlyTracker stored the calibration information.
% Default: 'calibration.mat'.
% processedtrxfilestr: Where to store the processed JAABA-compatible
% trajectories. Default: 'processed_trx.mat'. 
% processedflytrackerfilestr: Where to store the processed
% FlyTracker-compatible trajectories. Default: 'processed-movie-track.mat'.
% perframedirstr: Directory to store the wing features to. Default:
% 'perframe'. 
% indicatorfilestr: Where to store stimulation onset and offset times to.
% Default: indicatordata.mat. 
% metadatafilest: Where FlyBowlDataCapture stored experiment metadata.
% Default: Metadata.xml

function res = postProcessFlyTracker4JAABA(expdir,varargin)

codedir = fileparts(mfilename('fullpath'));
settingsdir = fullfile(codedir,'settings');

if ~exist('hmm_multiseq_1d','var'),
  addpath(fullfile(codedir,'hmm'));
end

maxFlyTrackerNanInterpFrames = 5;

[intrxfilestr,inflytrackerstr,...
  calibfilestr,proctrxfilestr,...
  procflytrackerfilestr,perframedirstr,...
  indicatorfilestr,...
  maxFlyTrackerNanInterpFrames,...
  metadatafilestr,...
  nmale,...
  leftovers] = ...
  myparse_nocheck(varargin,'trxfilestr','movie/movie_JAABA/trx.mat',...
  'flytrackerstr','movie/movie-track.mat',...
  'calibfilestr','calibration.mat',...
  'processedtrxfilestr','processed_trx.mat',...
  'processedflytrackerfilestr','processed-movie-track.mat',...
  'perframedirstr','perframe',...
  'indicatorfilestr','indicatordata.mat',...
  'maxFlyTrackerNanInterpFrames',maxFlyTrackerNanInterpFrames,...
  'metadatafilestr','Metadata.xml',...
  'nmale',[]);

intrxfile = fullfile(expdir,intrxfilestr);
assert(exist(intrxfile,'file')>0);
inflytrackerfile = fullfile(expdir,inflytrackerstr);
assert(exist(inflytrackerfile,'file')>0);
calibfile = fullfile(expdir,calibfilestr);
assert(exist(calibfile,'file')>0);
outtrxfile = fullfile(expdir,proctrxfilestr);
outflytrackerfile = fullfile(expdir,procflytrackerfilestr);
perframedir = fullfile(expdir,perframedirstr);
indicatorfile = fullfile(expdir,indicatorfilestr);

trxdataload = load(intrxfile);
trxin = trxdataload.trx;
cald = load(calibfile);
calib = cald.calib;
ftdataload = load(inflytrackerfile);
ftdatain = ftdataload.trk.data;

% split ids at nans that are more than maxFlyTrackerNanInterpFrames
% interpolate nans < maxFlyTrackerNanInterpFrames long
fns_check = {'x','y','theta','a','b'};
fns_perframe = {'x','y','theta','a','b',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','timestamps',...
  'xwingl','ywingl','xwingr','ywingr'};
anglefns = {'theta','theta_mm','ori',...
  'leg 1 ang','leg 2 ang','leg 3 ang','leg 4 ang','leg 5 ang','leg 6 ang',...
  'wing l ang','wing r ang'};

%% make row vectors and remove nans

trxout = [];
T = size(ftdatain,2);
d = size(ftdatain,3);
ftdataout = nan([0,T,d]);
sourceid = [];
ninterpframes = 0;
for i = 1:numel(trxin),
  trxcurr = trxin(i);
  % make row vectors
  for j = 1:numel(fns_perframe),
    fn = fns_perframe{j};
    trxcurr.(fn) = trxcurr.(fn)(:)';
  end

  ftdatacurr = ftdatain(i,:,:);
  isbad = false(1,trxcurr.nframes);
  for j = 1:numel(fns_check),
    isbad = isbad | isnan(trxcurr.(fns_check{j})(:)');
  end
  trxcurr.id = numel(trxout)+1;
  if ~any(isbad),
    trxout = structappend(trxout,trxcurr);
    sourceid(end+1) = i; %#ok<AGROW>
    ftdataout(trxcurr.id,:,:) = ftdatacurr;
    continue;
  end
  if all(isbad),
    continue;
  end
  % first frame bad, just crop
  if isbad(1),
    off = find(~isbad,1)-1;
    trxcurr.firstframe = trxcurr.firstframe+off;
    trxcurr.off = 1-trxcurr.firstframe;
    trxcurr.nframes = trxcurr.nframes-off;
    for j = 1:numel(fns_perframe),
      fn = fns_perframe{j};
      trxcurr.(fn) = trxcurr.(fn)(off+1:end);
    end
    isbad = isbad(off+1:end);
  end
  % last frame bad, just crop
  if isbad(end),
    off = trxcurr.nframes-find(~isbad,1,'last');
    trxcurr.endframe = trxcurr.endframe-off;
    trxcurr.nframes = trxcurr.nframes-off;
    for j = 1:numel(fns_perframe),
      fn = fns_perframe{j};
      trxcurr.(fn) = trxcurr.(fn)(1:end-off);
    end
    isbad = isbad(1:end-off);
  end
  if ~any(isbad),
    trxout = structappend(trxout,trxcurr);
    sourceid(end+1) = i; %#ok<AGROW>
    ftdataout(trxcurr.id,:,:) = ftdatacurr;
    continue;
  end
  % bad sequences in the middle
  [st,en] = get_interval_ends(isbad);
  badlength = en-st;
  idxinterp = find(badlength <= maxFlyTrackerNanInterpFrames);
  for jj = 1:numel(idxinterp),
    j = idxinterp(jj);
    t0 = st(j);
    t1 = en(j)-1;
    for k = 1:numel(fns_perframe),
      fn = fns_perframe{k};
      feat0 = trxcurr.(fn)(t0-1);
      feat1 = trxcurr.(fn)(t1+1);
      if ismember(fn,anglefns),
        featinterp = modrange(linspace(0,modrange(feat1-feat0,-pi,pi),badlength(j)+2)+feat0,-pi,pi);
      else
        featinterp = linspace(feat0,feat1,badlength(j)+2);
      end
      trxcurr.(fn)(t0-1:t1+1) = featinterp;
    end
    for k = 1:d,
      feat = ftdatacurr(:,:,k); 
      feat0 = feat(t0-1);
      feat1 = feat(t1+1);
      if ismember(fn,anglefns),
        featinterp = modrange(linspace(0,modrange(feat1-feat0,-pi,pi),badlength(j)+2)+feat0,-pi,pi);
      else
        featinterp = linspace(feat0,feat1,badlength(j)+2);
      end
      feat(t0-1:t1+1) = featinterp;
      ftdatacurr(:,:,k) = feat;
    end
    isbad(t0:t1) = false;
    ninterpframes = ninterpframes + (t1-t0+1);
  end
  isgood = ~isbad;
  [st,en] = get_interval_ends(isgood);
  for j = 1:numel(st),
    t0 = st(j);
    t1 = en(j)-1;
    trxnew = trxcurr;
    trxnew.firstframe = t0 + trxcurr.firstframe - 1;
    trxnew.endframe = t1 + trxcurr.firstframe - 1;
    trxnew.nframes = t1-t0+1;
    trxnew.off = 1-trxnew.firstframe;
    trxnew.id = numel(trxout)+1;
    for k = 1:numel(fns_perframe),
      fn = fns_perframe{k};
      trxnew.(fn) = trxcurr.(fn)(t0:t1);
    end
    trxout = structappend(trxout,trxnew);
    sourceid(end+1) = i; %#ok<AGROW>
    ftdatanew = nan(1,T,d);
    ftdatanew(:,t0:t1,:) = ftdatacurr(:,t0:t1,:);
    ftdataout(trxnew.id,:,:) = ftdatanew;
  end
  
end
  
%% center on arena center

arena = struct;
arena.x = calib.centroids(1) / calib.PPM;
arena.y = calib.centroids(2) / calib.PPM;
arena.r = calib.r / calib.PPM;

for i = 1:numel(trxout),
  trxout(i).x_mm = trxout(i).x_mm - arena.x; %#ok<AGROW> 
  trxout(i).y_mm = trxout(i).y_mm - arena.y; %#ok<AGROW> 
end

arena.x = 0;
arena.y = 0;

%% sex classification

metadatafile = fullfile(expdir,metadatafilestr);
metadata = ReadMetadataFile(metadatafile);

sex_classification_info = [];
if strcmpi(metadata.gender,'b'),
  if ~isempty(nmale),
    median_a_mm = nan(numel(trxout),1);
    for i = 1:numel(trxout),
      median_a_mm(i) = median(trxout(i).a_mm,'omitnan');
    end
    [~,sizeorder] = sort(median_a_mm);
    for i = 1:numel(trxout),
      if ismember(i,sizeorder(1:nmale)),
        trxout(i).sex = 'm'; %#ok<AGROW> 
      else
        trxout(i).sex = 'f'; %#ok<AGROW> 
      end
    end
  else
    [trx1,sex_classification_info] = FlyDiscoClassifySex(expdir,'trx',trxout,...
      'settingsdir',settingsdir,...
      'analysis_protocol','current_bubble',...
      'dosave',false,...
      'verbose',0);
    for i = 1:numel(trxout),
      nmalecurr = nnz(strcmp(trx1(i).sex,'M'));
      if nmalecurr >= trx1(i).nframes,
        trxout(i).sex = 'm'; %#ok<AGROW> 
      else
        trxout(i).sex = 'f'; %#ok<AGROW> 
      end
    end
  end
else
  for i = 1:numel(trxout),
    trxout(i).sex = metadata.gender; %#ok<AGROW>
  end

end

%% save output

trxdataload.trx = trxout;
save(outtrxfile,'-struct','trxdataload');
ftdataload.trk.data = ftdataout;
save(outflytrackerfile,'-struct','ftdataload');

fprintf('Postprocessed from %s->%s and\n%s->%s\nInterpolated %d frames\n%d trajectories -> %d trajectories. %d male, %d female.\n',...
  intrxfile,outtrxfile,inflytrackerfile,outflytrackerfile,ninterpframes,numel(trxin),numel(trxout),...
  nnz([trxout.sex]=='m'),nnz([trxout.sex]=='f'));

%% compute wing features

tdout = FlyTracker2WingTracking_helper(outflytrackerfile,outtrxfile,perframedir,outtrxfile,struct('fakectrax',false),arena);

res = struct;
res.trxout = tdout.trx;
res.ftdataout = ftdataout;
res.sex_classification_info = sex_classification_info;
res.processedtrxfile = outtrxfile;
res.processedflytrackerfile = outflytrackerfile;

%% determine stimulation onset and offsets

indicatorLED = estimateActivationTiming(expdir,leftovers{:});
if ~isempty(indicatorLED),
  save(indicatorfile,'indicatorLED');
end
res.indicatorLED = indicatorLED;

%% debug

if false,
  
  for i = 1:numel(trxout), %#ok<UNRCH> 
    
    id = sourceid(i);
    off1 = trxout(i).firstframe-trxin(id).firstframe;
    assert(all(isnan(trxin(id).x(off1+1:off1+trxout(i).nframes)) | trxout(i).x == trxin(id).x(off1+1:off1+trxout(i).nframes)));
    assert(numel(trxout(i).x)==trxout(i).nframes);
    assert((trxout(i).endframe-trxout(i).firstframe+1)==trxout(i).nframes);
    assert(all(all(isnan(ftdataout(i,:,:)) | ftdataout(i,:,:) == ftdatain(id,:,:))));
    
  end
  
  clf;
  hax(1) = subplot(6,1,1);
  hold on;
  for i = 1:numel(trxin),
    plot(trxin(i).firstframe:trxin(i).endframe,trxin(i).y,'.-');
  end
  hax(2) = subplot(6,1,2);
  hold on;
  for i = 1:numel(trxout),
    plot(trxout(i).firstframe:trxout(i).endframe,trxout(i).y,'.-');
  end
  hax(3) = subplot(6,1,3);
  hold on;
  for i = 1:numel(trxout),
    plot(trxout(i).firstframe:trxout(i).endframe,trxout(i).theta,'.-');
  end
  hax(4) = subplot(6,1,4);
  hold on;
  j = find(strcmp(ftdataload.trk.names,'pos x'));
  for i = 1:size(ftdataout,1),
    plot(ftdataout(i,:,j),'.-');
  end
  hax(5) = subplot(6,1,5);
  hold on;
  j = find(strcmp(ftdataload.trk.names,'pos y'));
  for i = 1:size(ftdataout,1),
    plot(ftdataout(i,:,j),'.-');
  end
  hax(6) = subplot(6,1,6);
  hold on;
  j = find(strcmp(ftdataload.trk.names,'ori'));
  for i = 1:size(ftdataout,1),
    plot(ftdataout(i,:,j),'.-');
  end
  linkaxes(hax,'x');
end