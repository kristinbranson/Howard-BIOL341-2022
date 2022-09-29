function [trx,summary_diagnostics,X] = FlyDiscoClassifySex(expdir,varargin)

version = '0.2';

%% parse parameters
[analysis_protocol,...
  settingsdir,...
  datalocparamsfilestr,...
  outdir,...
  dosave,...
  override_gender,...
  trxfile,...
  trx,...
  verbose] = myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'outdir',[],...
  'dosave',true,...
  'override_gender','',...
  'trxfile','',...
  'trx',[],...
  'verbose',1);

if ~isempty(trx) && ~isempty(trxfile),
  fprintf('Classifying sex for %s\n',expdir);
end

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% %%
% logger = PipelineLogger(expdir,mfilename(),...
%         dataloc_params,'classifysex_logfilestr',...
%         settingsdir,analysis_protocol,'versionstr',version);     

%% load the data
fprintf('Loading data...\n');

if isempty('trx'),
  if isempty(trxfile),
    trxfile = fullfile(expdir,dataloc_params.trxfilestr);
  end
  load(trxfile,'trx');
end
% first and end frame
firstframe = min([trx.firstframe]); 
endframe = max([trx.endframe]);

%% data locations
sexclassifierparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.sexclassifierparamsfilestr);
sexclassifierintxtfile = fullfile(settingsdir,analysis_protocol,dataloc_params.sexclassifiertxtfilestr);
sexclassifier_params = ReadParams(sexclassifierparamsfile);
sexclassifierin = ReadParams(sexclassifierintxtfile);

%% output files
if dosave,
  sexclassifieroutmatfile = fullfile(expdir,outdir,...
    dataloc_params.sexclassifiermatfilestr);
  sexclassifierdiagnosticsfile = fullfile(expdir,outdir,...
    dataloc_params.sexclassifierdiagnosticsfilestr);
  trxfileout = fullfile(expdir,outdir,dataloc_params.trxfilestr);
  tftrxappend = strcmp(trxfile,trxfileout);
end

%% compute areas
fprintf('Computing area...\n');

nflies = numel(trx); 

if nflies <= 0,
  error('No flies tracked in this video');
end

X = cell(1,nflies);
for fly = 1:nflies,
  % compute area
  area = (2*trx(fly).a_mm).*(2*trx(fly).b_mm)*pi;
  badidx = isinf(area) | isnan(area);
  if any(isnan(area)),
    warning('FlyBubbleClasisfySex:area',...
      'NaNs found in area for fly %d of experiment %s',fly,expdir);
    area(badidx) = inf;
  end
  areacurr = SmoothAreaOutliers(area,...
    sexclassifierin.areasmooth_filterorder,...
    sexclassifierin.areasmooth_maxfreq,...
    sexclassifierin.areasmooth_maxerrx);
  X{fly} = areacurr(:);
end

%% fix sex if gender is not 'b'

% read gender
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
metadata = ReadMetadataFile(metadatafile);
if ~isempty(override_gender),
  metadata.gender = override_gender;
end

% Fix possible metadata error
if isfield(metadata, 'gender') && strcmpi(metadata.gender, 'both') ,
  metadata.gender = 'b' ;
end

if isfield(metadata,'gender') && ~strcmpi(metadata.gender,'b'),
  fprintf('gender is not "b", not doing sex classification, just setting sex to %s for all flies',upper(metadata.gender));

  % set sex to metadata.gender for all flies
  % also set diagnostics that we can
  clear diagnostics
  mean_area_all = nanmean(cat(1,X{:}));
  var_area_all = nanstd(cat(1,X{:}),1);
  for fly = 1:numel(trx),
    trx(fly).sex = repmat({upper(metadata.gender)},[1,trx(fly).nframes]);
    diagnostics_curr = struct;
    diagnostics_curr.normhmmscore = nan;
    diagnostics_curr.nswaps = 0;
    diagnostics_curr.meanabsdev = nanmean(abs(X{fly}-mean_area_all));
    diagnostics(fly) = diagnostics_curr; %#ok<AGROW>
  end
  
  mu_area = nan(1,2);
  var_area = nan(1,2);
  if strcmpi(metadata.gender,'f'),
    mu_area(2) = mean_area_all;
    var_area(2) = var_area_all;
  elseif strcmpi(metadata.gender,'m'),
    mu_area(1) = mean_area_all;
    var_area(1) = var_area_all;
  end
  state2sex = {'M','F'};
  ll = [];
  
else
  
  %% learn a 2-state HMM for area in an unsupervised manner
  
  fprintf('gender = "b", learning 2-state HMM...\n');
  
  % initialize parameters
  nstates = 2;
  if isfield(sexclassifierin,'ptrans')
    fprintf(1,'''ptrans'' field present in SC settings...\n');

    ptrans = sexclassifierin.ptrans;
    psame = 1-ptrans;
  else
    fprintf(1,'''ptrans'' field not present in SC settings. Using ''psame'' field...\n');
   
    psame = sexclassifierin.psame;
    ptrans = 1-psame;
  end
  fprintf(1,' ... setting [psame ptrans] = [%.03g %.03g]\n',psame,ptrans);
  ptransmat = ones(nstates)*ptrans;
  ptransmat(eye(nstates)==1) = psame;
  
  if isfield(sexclassifier_params,'frac_female'),
    prior = [1-sexclassifier_params.frac_female,sexclassifier_params.frac_female];
  else
    prior = ones(1,nstates)/nstates;
  end
  state2sex = cell(1,nstates);
  
  % em for hmm
  [mu_area,var_area,ll,mu_area_km,var_area_km]=hmm_multiseq_1d(X,nstates,ptransmat,...
    sexclassifier_params.niters_em,sexclassifier_params.tol_em,...
    [],prior,verbose);
  if mu_area(1) > mu_area(2),
    mu_area = mu_area(end:-1:1);
    var_area = var_area(end:-1:1);
  end
  [~,maxi] = max(mu_area);
  state2sex{maxi} = 'F';
  [~,mini] = min(mu_area);
  state2sex{mini} = 'M';
  
  if verbose > 0,
    fprintf(1,'mu_area(1) mu_area_km(1) mu_area(2) mu_area_km(2): %0.3f %0.3f %0.3f %0.3f\n',...
      mu_area(1),mu_area_km(1),mu_area(2),mu_area_km(2));
    fprintf(1,'vr_area(1) vr_area_km(1) vr_area(2) vr_area_km(2): %0.3f %0.3f %0.3f %0.3f\n',...
      var_area(1),var_area_km(1),var_area(2),var_area_km(2));
  end    
    
  %% save classifier
  
  filterorder = sexclassifierin.areasmooth_filterorder; 
  maxfreq = sexclassifierin.areasmooth_maxfreq; 
  maxerrx = sexclassifierin.areasmooth_maxerrx; 
  
  if dosave,
    try
      if exist(sexclassifieroutmatfile,'file'),
        delete(sexclassifieroutmatfile);
      end
    catch ,
      % Ignore any errors that occur
    end
    try
      save(sexclassifieroutmatfile,'mu_area','var_area','ptrans','prior','ll',...
        'nstates','state2sex','maxerrx','maxfreq','filterorder','version','analysis_protocol');
    catch ME,
      warning('FlyDiscoClassifySex:save',...
        'Could not save to file %s: %s',sexclassifieroutmatfile,getReport(ME));
      fprintf('Could not save to file %s: %s\n',sexclassifieroutmatfile,getReport(ME));
    end
  end
  
  %% classify sex
  
  clear diagnostics;
  for fly = 1:numel(trx),
    
    % Viterbi to classify per-frame
    [trx(fly).sex,diagnostics(fly)] = ClassifySex(X{fly}',...
      mu_area,var_area,ptransmat,state2sex); %#ok<AGROW>
    
  end
  
end

%% count number of flies, females, males
counts = struct;
counts.nfemales = zeros(1,endframe);
counts.nmales = zeros(1,endframe);
counts.nflies = zeros(1,endframe);
for fly = 1:numel(trx),
  isfemale = strcmp(trx(fly).sex,'F');
  ismale = strcmp(trx(fly).sex,'M');
  counts.nfemales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nfemales(trx(fly).firstframe:trx(fly).endframe) + double(isfemale);
  counts.nmales(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nmales(trx(fly).firstframe:trx(fly).endframe) + double(ismale);
  counts.nflies(trx(fly).firstframe:trx(fly).endframe) = ...
    counts.nflies(trx(fly).firstframe:trx(fly).endframe) + 1;
end

% ignore part of video with no flies
fns = fieldnames(counts);
for i = 1:numel(fns),
  fn = fns{i};
  counts.(fn)(1:firstframe-1) = nan;
  counts.(fn)(endframe+1:end) = nan;
end

%% summary diagnostics

summary_diagnostics = struct;
ifemale = find(strcmp(state2sex,'F'));
imale = find(strcmp(state2sex,'M'));
summary_diagnostics.classifier_mu_area_female = mu_area(ifemale);
summary_diagnostics.classifier_mu_area_male = mu_area(imale);
summary_diagnostics.classifier_mu_area_female = mu_area(ifemale);
summary_diagnostics.classifier_mu_area_male = mu_area(imale);
summary_diagnostics.classifier_var_area_female = var_area(ifemale);
summary_diagnostics.classifier_var_area_male = var_area(imale);
if isempty(ll),
  summary_diagnostics.classifier_loglik = nan;
else
  summary_diagnostics.classifier_loglik = ll(end);
end
summary_diagnostics.classifier_niters = numel(ll);

fns = fieldnames(diagnostics);
for i = 1:numel(fns),
  fn = fns{i};
  summary_diagnostics.(sprintf('mean_%s',fn)) = nanmean([diagnostics.(fn)]);
  summary_diagnostics.(sprintf('median_%s',fn)) = nanmedian([diagnostics.(fn)]);
  summary_diagnostics.(sprintf('std_%s',fn)) = nanstd([diagnostics.(fn)],1);
  summary_diagnostics.(sprintf('min_%s',fn)) = min([diagnostics.(fn)]);
  summary_diagnostics.(sprintf('max_%s',fn)) = max([diagnostics.(fn)]);
end

fns = fieldnames(counts);
for i = 1:numel(fns),
  fn = fns{i};
  summary_diagnostics.(sprintf('mean_%s',fn)) = nanmean([counts.(fn)]);
  summary_diagnostics.(sprintf('median_%s',fn)) = nanmedian([counts.(fn)]);
  summary_diagnostics.(sprintf('std_%s',fn)) = nanstd([counts.(fn)],1);
  summary_diagnostics.(sprintf('min_%s',fn)) = min([counts.(fn)]);
  summary_diagnostics.(sprintf('max_%s',fn)) = max([counts.(fn)]);
end

fns1 = fieldnames(summary_diagnostics);
if verbose > 0,
  fprintf('Summary diagnostics:\n');
  for i = 1:numel(fns1),
    fprintf('  %s: %f\n',fns1{i},summary_diagnostics.(fns1{i}));
  end
end

%% write diagnostics

% sexclassifierinfo = logger.runInfo;
% sexclassifierinfo.version = version;

if dosave,
  if exist(sexclassifierdiagnosticsfile,'file'),
    try %#ok<TRYNC>
      delete(sexclassifierdiagnosticsfile);
    end
  end      
  fid = fopen(sexclassifierdiagnosticsfile,'w');
  if fid < 0,
    warning('FlyDiscoClassifySex:diags',...
      'Could not open file %s for writing, printing diagnostics to stdout instead',sexclassifierdiagnosticsfile);
    fid = 1;
  end
else
  fid = 1;
end
if dosave || verbose > 0,
  fns = fieldnames(summary_diagnostics);
  for i = 1:numel(fns),
    fprintf(fid,'%s,%f\n',fns{i},summary_diagnostics.(fns{i}));
  end
end
% fns = fieldnames(sexclassifierinfo);
% for i = 1:numel(fns),
%   val = sexclassifierinfo.(fns{i});
%   if ischar(val)
%     fprintf(fid,'%s,%s\n',fns{i},val);
%   end
% end

if dosave && fid > 1,
  fclose(fid);
end

%% resave

if dosave,
  try
    if tftrxappend
      save('-append',trxfileout,'trx');
    else
      tmp = load(trxfile);
      tmp.trx = trx;
      %tmp.sexclassifierinfo = sexclassifierinfo;
      save(trxfileout,'-struct','tmp');
    end
  catch
    try
      tmp = load(trxfile);
      trxfilebak = [trxfileout,'.bak'];
      unix(sprintf('mv %s %s',trxfileout,trxfilebak));
      tmp.trx = trx;
%       tmp.sexclassifierinfo = sexclassifierinfo;
      save(trxfileout,'-struct','tmp');
    catch ME,
      warning('FlyDiscoClassifySex:save',...
        'Could not save to file %s: %s',trxfileout,getReport(ME));
      fprintf('Could not save to file %s: %s',trxfileout,getReport(ME));
    end
  end
end

%logger.close();
