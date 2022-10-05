% data = LoadTracking(expdirs,...)
%
% Load tracking information from experiment directories expdirs. 
% 
% Input:
% expdirs: cell of full paths to experiment directories. 
% 
% Output: 
% data: data is a struct. It has a field for each of the following pieces
% of trajectory information:
% x_mm: x-coordinate of centroid in millimeters, ntrajs x T array.
% y_mm: y-coordinate of centroid in millimeters, ntrajs x T array.
% a_mm: quarter-major axis length in millimeters, ntrajs x T array.
% b_mm: quarter-minor axis length in millimeters, ntrajs x T array.
% x_px: x-coordinate of centroid in pixels, ntrajs x T array.
% y_px: y-coordinate of centroid in pixels, ntrajs x T array.
% theta_rad: orientation of fly in radians, ntrajs x T array.
% a_px: quarter-major axis length in pixels, ntrajs x T array.
% b_px: quarter-minor axis length in pixels, ntrajs x T array.
% xwingl_px: x-coordinate of left wing tip, ntrajs x T array.
% ywingl_px: y-coordinate of left wing tip, ntrajs x T array.
% xwingr_px: x-coordinate of right wing tip, ntrajs x T array.
% ywingr_px: y-coordinate of right wing tip, ntrajs x T array.
% wing_anglel_rad: left wing angle in radians, ntrajs x T array. 
% wing_angler_rad: right wing angle in radians, ntrajs x T array. 
% timestamp: when the frame was recorded (same for all flies in a video),
% ntrajs x T array.
% Here, ntrajs is the total number of fly trajectories in all experiment
% directories. T is the maximum number of frames in any video. 
% data also has a field metadata which is a table with the following
% information: 
% exp_name: Name of experiment (same for all flies in a video), string. 
% exp_num: Index for the experiment (same for all flies in a video), number
% 1, 2,... 
% pxpermm: Pixels per millimeter conversion factor (same for all flies in
% a video), number.
% fps: Frames per second conversion factor (same for all flies in
% a video), number. 
% id: Index for the fly in the experiment, number 1, 2, ..., 10
% sex: 'f' or 'm' indicating whether the fly is female or male. 
% activation_startframes: onsets of activation periods (same for all flies
% in a video), 1 x nactivation array. 
% activation_endframes: offsets of activation periods (same for all flies
% in a video), 1 x nactivation array.
% activation_intensities: intensity of red light during activation periods
% (same for all flies in a video), 1 x nactivation array.
% activation_pulsewidths: width of pulses of light during activation
% periods in us (same for all flies in a video), 1 x nactivation array.
% activation_pulseperiods: period for pulsing light during activation
% periods in us (same for all flies in a video), 1 x nactivation array. 
% Here, nactivation is the number of activation periods for the video. 
%
%
% Example ways to access data:
% data = LoadTracking(expdirs)
% x = data.x_mm(data.metadata.exp_num==1,:);
% Returns a matrix of size ntraj x T
% y = data.y_mm(17:18,101:110)
% Returns a matrix of size 2 x 10

function data = LoadTracking(expdirs,varargin)

[trxfilestr,...
  indicatorfilestr] = ...
  myparse(varargin,...
  'trxfilestr','processed_trx.mat',...
  'indicatorfilestr','indicatordata.mat');

if ischar(expdirs),
  expdirs = {expdirs};
end

nexps = numel(expdirs);

copy_perframe = {'x_mm','y_mm','a_mm','b_mm'};
outnames_perframe =  {'x_px','y_px','theta_rad','a_px','b_px','xwingl_px','ywingl_px','xwingr_px','ywingr_px','wing_anglel_rad','wing_angler_rad'};
innames_perframe = {'x',   'y',   'theta',    'a',   'b',   'xwingl',   'ywingl',   'xwingr',   'ywingr'   ,'wing_anglel',    'wing_angler'};
outnames_perframe = [copy_perframe,outnames_perframe];
innames_perframe = [copy_perframe,innames_perframe];

sexes = {'f','m'};

data.summary.exps = [];
data.summary.flies = [];
data.summary.activation = [];
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

  expcurr = struct;
  expinfo = struct;
  expinfo.expname = string(expname);
  expinfo1 = parseExpName(expname);
  fns = fieldnames(expinfo1);
  for j = 1:numel(fns),
    expinfo.(fns{j}) = string(expinfo1.(fns{j}));
  end
  expinfo.expname = string(expname);
  expinfo.pxpermm = td.trx(1).pxpermm;
  expinfo.fps = td.trx(1).fps;
  expinfo.nframes = max([td.trx.endframe]);
  expinfo.nflies = numel(td.trx);
  expcurr.summary = struct2table(expinfo);
  if expi == 1,
    data.summary.exps = repmat(expcurr.summary,[nexps,1]);
  else
    data.summary.exps(expi,:) = expcurr.summary;
  end

  for flyi = 1:numel(td.trx),

    flycurr = struct;
    flycurr.expnum = expi;
    flycurr.flynum = flyi;
    flycurr.sex = categorical({td.trx(flyi).sex},sexes);
    flycurr.startframe = td.trx(flyi).firstframe;
    flycurr.endframe = td.trx(flyi).endframe;
    flycurr.nframes = td.trx(flyi).nframes;
    tabcurr = struct2table(flycurr);
    data.summary.flies = [data.summary.flies;tabcurr];
    flycurr = rmfield(flycurr,{'expnum','flynum'});

    % tracking
    t0 = td.trx(flyi).firstframe;
    t1 = td.trx(flyi).endframe;
    for inti = 1:numel(outnames_perframe),
      flycurr.(outnames_perframe{inti})(t0:t1) = td.trx(flyi).(innames_perframe{inti});
    end

    if flyi == 1,
      expcurr.fly = repmat(flycurr,[numel(td.trx),1]);
    else
      expcurr.fly(flyi) = flycurr;
    end

  end

  % activation
  expcurr.activation = [];
  if isopto,
    for inti = 1:numel(ind.indicatorLED.startframe),
      row = [expi, inti, ind.indicatorLED.startframe(inti), ind.indicatorLED.endframe(inti), ...
        ind.indicatorLED.intensity(inti), ind.indicatorLED.pulsewidths(inti), ...
        ind.indicatorLED.pulseperiods(inti)];
      tab = array2table(row,'VariableNames',{'expnum','intervalnum','startframe','endframe','intensity','pulsewidth','pulseperiod'});
      data.summary.activation = [data.summary.activation; tab];
      tab = array2table(row(3:end),'VariableNames',{'startframe','endframe','intensity','pulsewidth','pulseperiod'});
      expcurr.activation = [expcurr.activation; tab];
    end
  end

  expcurr.timestamps = td.timestamps;
  if expi == 1,
    data.exp = repmat(expcurr,[nexps,1]);
  else
    data.exp(expi) = expcurr;
  end


end
    