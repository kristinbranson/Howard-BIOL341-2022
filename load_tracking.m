function tab = load_tracking(expdirs,varargin)

[trxfilestr,...
  indicatorfilestr] = ...
  myparse(varargin,...
  'trxfilestr','processed_trx.mat',...
  'indicatorfilestr','indicatordata.mat');

if ischar(expdirs),
  expdirs = {expdirs};
end

nexps = numel(expdirs);

tablenames_perexp = {'exp_name','exp_num'};
tablenames_perframe =  {'x_px','y_px','theta_rad','a_px','b_px','xwingl_px','ywingl_px','xwingr_px','ywingr_px','wing_anglel_rad','wing_angler_rad','timestamp'};
structnames_perframe = {'x',   'y',   'theta',    'a',   'b',   'xwingl',   'ywingl',   'xwingr',   'ywingr'   ,'wing_anglel',    'wing_angler',    'timestamps'};
copy_perframe = {'x_mm','y_mm','a_mm','b_mm'};
tablenames_perframe = [tablenames_perframe,copy_perframe];
structnames_perframe = [structnames_perframe,copy_perframe];
tablenames_perfly = {'id','sex'};
structnames_activation = {'startframe','endframe','intensity','pulsewidths','pulseperiods'};
tablenames_activation = {'activation_startframes','activation_endframes','activation_intensities','activation_pulsewidths','activation_pulseperiods'};

tablenames = [tablenames_perexp,tablenames_perfly,tablenames_perframe,{'frame'},tablenames_activation];

sexes = {'f','m'};

maxT = 0;
ntrajs = nan(nexps,1);
for expi = 1:nexps,
  expdir = expdirs{expi};
  trxfile = fullfile(expdir,trxfilestr);
  td = matfile(trxfile);
  T = numel(td.timestamps);
  maxT = max(maxT,T);
  ntrajs(expi) = numel(td.trx);
end

rowidx = 1;
tab = [];
for expi = 1:nexps,
  expdir = expdirs{expi};
  [~,expname] = fileparts(expdir);
  trxfile = fullfile(expdir,trxfilestr);
  td = load(trxfile);

  indicatorfile = fullfile(expdir,indicatorfilestr);
  isopto = false;
  if exist(indicatorfile,'file'),
    ind = load(indicatorfile);
    isopto = isfield(ind,'indicatorLED') && ~isempty(ind.indicatorLED);
  end

  for i = 1:numel(td.trx),
    row = cell(1,numel(tablenames));
    idx = 1;

    % per-experiment
    row{idx} = expname;
    idx = idx + 1;
    row{idx} = expi;
    idx = idx + 1;
    % per-fly
    row{idx} = i;
    idx = idx + 1;
    row{idx} = categorical({td.trx(i).sex},sexes);
    idx = idx + 1;

    % tracking
    t0 = td.trx(i).firstframe;
    t1 = td.trx(i).endframe;
    for j = 1:numel(tablenames_perframe),
      x = nan(1,maxT);
      x(t0:t1) = td.trx(i).(structnames_perframe{j});
      row{idx} = x;
      idx = idx + 1;
    end
    row{idx} = 1:maxT;
    idx = idx + 1;

    % activation
    for j = 1:numel(tablenames_activation),
      if isopto,
        row{idx} = ind.indicatorLED.(structnames_activation{j});
      else
        row{idx} = [];
      end
      idx = idx + 1;
    end

    tabcurr = cell2table(row,'VariableNames',tablenames);
    if isempty(tab),
      tab = repmat(tabcurr,[sum(ntrajs),1]);
    else
      tab(rowidx,:) = tabcurr; %#ok<AGROW> 
    end
    rowidx = rowidx + 1;
  end

end
    