function res = postProcessFlyTracker4JAABA(expdir,varargin)

maxFlyTrackerNanInterpFrames = 5;
thresh_female_a_mm = .64; % min quarter major axis length for female flies

[intrxfilestr,inflytrackerstr,...
  calibfilestr,proctrxfilestr,...
  procflytrackerfilestr,perframedirstr,...
  indicatorfilestr,...
  maxFlyTrackerNanInterpFrames,...
  thresh_female_a_mm,metadatafilestr,...
  leftovers] = ...
  myparse_nocheck(varargin,'trxfilestr','movie/movie_JAABA/trx.mat',...
  'flytrackerstr','movie/movie-track.mat',...
  'calibfilestr','calibration.mat',...
  'processedtrxfilestr','processed_trx.mat',...
  'processedflytrackerfilestr','processed-movie-track.mat',...
  'perframedirstr','perframe',...
  'indicatorfilestr','indicatordata.mat',...
  'maxFlyTrackerNanInterpFrames',maxFlyTrackerNanInterpFrames,...
  'thresh_female_a_mm',thresh_female_a_mm,...
  'metadatafilestr','Metadata.xml');

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


trxout = [];
T = size(ftdatain,2);
d = size(ftdatain,3);
ftdataout = nan([0,T,d]);
sourceid = [];
ninterpframes = 0;
for i = 1:numel(trxin),
  trxcurr = trxin(i);
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

metadatafile = fullfile(expdir,metadatafilestr);
metadata = ReadMetadataFile(metadatafile);

for i = 1:numel(trxout),
  trxout(i).dt = diff(trxout(i).timestamps); %#ok<AGROW>
  % center on arena center
  trxout(i).x_mm = trxout(i).x_mm - calib.centroids(1)/calib.PPM; %#ok<AGROW>
  trxout(i).y_mm = trxout(i).y_mm - calib.centroids(2)/calib.PPM; %#ok<AGROW>
  for j = 1:numel(fns_perframe),
    fn = fns_perframe{j};
    trxout(i).(fn) = trxout(i).(fn)(:)';
  end
  if strcmpi(metadata.gender,'b'),
    median_a_mm = median(trxout(i).a_mm,'omitnan');
    if median_a_mm >= thresh_female_a_mm,
      trxout(i).sex = 'f'; %#ok<AGROW>
    else
      trxout(i).sex = 'm'; %#ok<AGROW>
    end
  else
    trxout(i).sex = metadata.gender; %#ok<AGROW> 
  end
end

trxdataload.trx = trxout;
save(outtrxfile,'-struct','trxdataload');
ftdataload.trk.data = ftdataout;
save(outflytrackerfile,'-struct','ftdataload');

fprintf('Postprocessed from %s->%s and\n%s->%s\nInterpolated %d frames\n%d trajectories -> %d trajectories. %d male, %d female.\n',...
  intrxfile,outtrxfile,inflytrackerfile,outflytrackerfile,ninterpframes,numel(trxin),numel(trxout),...
  nnz([trxout.sex]=='m'),nnz([trxout.sex]=='f'));

arena = struct;
arena.x = calib.centroids(1);
arena.y = calib.centroids(2);
arena.r = calib.r;

tdout = FlyTracker2WingTracking_helper(outflytrackerfile,outtrxfile,perframedir,outtrxfile,struct('fakectrax',false),arena);

res = struct;
res.trxout = tdout.trx;
res.ftdataout = ftdataout;
res.processedtrxfile = outtrxfile;
res.processedflytrackerfile = outflytrackerfile;

indicatorLED = estimateActivationTiming(expdir,leftovers{:});
save(indicatorfile,'indicatorLED');
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