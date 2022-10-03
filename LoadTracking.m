% exptable = LoadTracking(expdirs,...)
%
% Load tracking information from experiment directories expdirs into a
% table. The table will have a row for each fly trajectory. It will have a
% column for each of the following variables:
% exp_name: Name of experiment (same for all flies in an experiment
% directory), string. 
% exp_num: Index for the experiment (same for all flies in an experiment
% directory), number 1, 2,... 
% pxpermm: Pixels per millimeter conversion factor (same for all flies in
% a video), number.
% fps: Frames per second conversion factor (same for all flies in
% a video), number. 
% id: Index for the fly in the experiment, number 1, 2, ..., 10
% sex: 'f' or 'm' indicating whether the fly is female or male. 
% x_mm: x-coordinate of centroid in millimeters, 1 x T array.
% y_mm: y-coordinate of centroid in millimeters, 1 x T array.
% a_mm: quarter-major axis length in millimeters, 1 x T array.
% b_mm: quarter-minor axis length in millimeters, 1 x T array.
% x_px: x-coordinate of centroid in pixels, 1 x T array.
% y_px: y-coordinate of centroid in pixels, 1 x T array.
% theta_rad: orientation of fly in radians, 1 x T array.
% a_px: quarter-major axis length in pixels, 1 x T array.
% b_px: quarter-minor axis length in pixels, 1 x T array.
% xwingl_px: x-coordinate of left wing tip, 1 x T array.
% ywingl_px: y-coordinate of left wing tip, 1 x T array.
% xwingr_px: x-coordinate of right wing tip, 1 x T array.
% ywingr_px: y-coordinate of right wing tip, 1 x T array.
% wing_anglel_rad: left wing angle in radians, 1 x T array. 
% wing_angler_rad: right wing angle in radians, 1 x T array. 
% timestamp: when the frame was recorded (same for all flies), 1 x T array.
% activation_startframes: onsets of activation periods (same for all flies
% in a video), 1 x nactivation. 
% activation_endframes: offsets of activation periods (same for all flies
% in a video), 1 x nactivation. 
% activation_intensities: intensity of red light during activation periods
% (same for all flies in a video), 1 x nactivation. 
% activation_pulsewidths: width of pulses of light during activation
% periods in us (same for all flies in a video), 1 x nactivation. 
% activation_pulseperiods: period for pulsing light during activation
% periods in us (same for all flies in a video), 1 x nactivation. 
% Above, T is the maximum number of frames in any video. 
% nactivation is the number of activation periods. 
%
% Input:
% expdirs: cell of full paths to experiment directories. 
%
% Example ways to access data:
% exptable = LoadTracking(expdirs)
% x = exptable.x_mm(exptable.exp_num==1,:);
% Returns a matrix of size ntraj x T
% y = exptable.y_mm(17:18,101:110)
% Returns a matrix of size 2 x 10

function tab = LoadTracking(expdirs,varargin)

[trxfilestr,...
  indicatorfilestr] = ...
  myparse(varargin,...
  'trxfilestr','processed_trx.mat',...
  'indicatorfilestr','indicatordata.mat');

if ischar(expdirs),
  expdirs = {expdirs};
end

nexps = numel(expdirs);

tablenames_perexp = {'exp_name','exp_num','pxpermm','fps'};
tablenames_perframe =  {'x_px','y_px','theta_rad','a_px','b_px','xwingl_px','ywingl_px','xwingr_px','ywingr_px','wing_anglel_rad','wing_angler_rad','timestamp'};
structnames_perframe = {'x',   'y',   'theta',    'a',   'b',   'xwingl',   'ywingl',   'xwingr',   'ywingr'   ,'wing_anglel',    'wing_angler',    'timestamps'};
copy_perframe = {'x_mm','y_mm','a_mm','b_mm'};
tablenames_perframe = [copy_perframe,tablenames_perframe];
structnames_perframe = [copy_perframe,structnames_perframe];
tablenames_perfly = {'id','sex'};
structnames_activation = {'startframe','endframe','intensity','pulsewidths','pulseperiods'};
tablenames_activation = {'activation_startframes','activation_endframes','activation_intensities','activation_pulsewidths','activation_pulseperiods'};

tablenames = [tablenames_perexp,tablenames_perfly,tablenames_perframe,tablenames_activation];

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
    row{idx} = td.trx(1).pxpermm;
    idx = idx + 1;
    row{idx} = td.trx(1).fps;
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
    