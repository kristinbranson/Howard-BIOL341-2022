% data = LoadTracking(expdirs,...)
%
% Load tracking information from experiment directories expdirs. 
% 
% Input:
% expdirs: cell of full paths to experiment directories. 
% 
% Output: 
% data: data is a struct. 
%
% data.summary contains the following summary information about all the
% videos loaded in:
%
% data.summary.exps is a table with a row for each video, describing all
% the experiments represented by this struct. It contains the following:
% data.summary.exps.expname(expi): Name of experiment, parsed from file
% path (string).
% data.summary.exps.exptype(expi): Type of experiment, parsed from expname.
% data.summary.exps.condition(expi): String indicating the flies'
% condition, parsed from experiment name (string). 
% data.summary.exps.rig(expi): String indicating which rig the experiment
% was recorded in (string). 
% data.summary.exps.timestamp(expi): String indicating when the experiment
% was recorded (string).
% data.summary.exps.pxpermm(expi): Pixels per millimeter conversion factor
% (number). 
% data.summary.exps.nframes(expi): Last frame tracked for any fly in the
% video (number). 
% data.summary.exps.nflies(expi): Number of fly trajectories for this video
% (number). 
%
% data.summary.flies is a table with a row for each fly in any video,
% describing the flies represented by this struct. It contains the
% following:
% data.summary.flies.expnum(flyi): Which experiment this fly is part of
% (number). 
% data.summary.flies.flynum(flyi): Which fly in experiment expnum this row
% corresponds to (number). 
% data.summary.flies.sex(flyi):  'f' or 'm' indicating whether the fly is
% female or male (categorical). 
% data.summary.flies.startframe(flyi): The first frame of the video this
% fly is tracked (number). 
% data.summary.flies.endframe(flyi): The last frame of the video this
% fly is tracked (number). 
% data.summary.flies.nframes(flyi): The number of video frames this
% fly is tracked for (number). 
% 
% data.summary.activation is a table with a row for each activation period
% in any video, describing the activation periods represented by this
% struct. It contains the following:
% data.summary.activation.expnum(inti): Which experiment this activation
% period is a part of (number).
% data.summary.activation.intervalnum(inti): Which activation interval in
% the video this corresponds to (number).
% data.summary.activation.startframe(inti): Which frame of the video this
% activation interval starts on (number). 
% data.summary.activation.endframe(inti): Which frame of the video this
% activation interval ends on (number). 
% data.summary.activation.intensity(inti): Intensity (%) of red light during
% activation periods (number).
% data.summary.activation.pulsewidth(inti): Width of pulses of light during
% activation, in microseconds (number). 
% data.summary.activation.pulseperiod(inti): Period of pulses of light during
% activation, in microseconds (number). 
% 
% data.exp is an array of structs, where data.exp(expnum) is a struct
% containing information related to experiment expnum. 
% data.exp(expnum).summary is a table with a single row containing
% information about this experiment. It is the same as row expnum of
% data.summary.exps. 
% data.exp(expnum).summary.expname: Name of experiment, parsed from file
% path (string).
% data.exp(expnum).summary.exptype: Type of experiment, parsed from expname.
% data.exp(expnum).summary.condition: String indicating the flies'
% condition, parsed from experiment name (string). 
% data.exp(expnum).summary.rig: String indicating which rig the experiment
% was recorded in (string). 
% data.exp(expnum).summarytimestamp: String indicating when the experiment
% was recorded (string).
% data.exp(expnum).summary.pxpermm: Pixels per millimeter conversion factor
% (number). 
% data.exp(expnum).summary.nframes: Last frame tracked for any fly in the
% video (number). 
% data.exp(expnum).summary.nflies: Number of fly trajectories for this video
% (number). 
%
% data.exp(expnum).fly is an array of structs, were
% data.exp(expnum).fly(flynum) contains trajectory data for fly flynum in
% experiment expnum. It contains the following fields:
%
% data.exp(expnum).fly(flynum).sex: 'f' or 'm' indicating whether the fly is
% female or male (categorical). 
% data.exp(expnum).fly(flynum).startframe: The first frame of the video
% this fly is tracked (number).
% data.exp(expnum).fly(flynum).endframe: The last frame of the video this
% fly is tracked (number). 
% data.exp(expnum).fly(flynum).nframes: The number of video frames this
% fly is tracked for (number). 
% The following fields are arrays, with an entry for each of the T =
% data.exp(expnum).nframes frames in the experiment.
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
%
% Example ways to access data:
%
% data = LoadTracking(expdirs)
% 
% Access data for a specific experiment and fly:
% x = data.exp(expnum).fly(flynum).x_mm
% Returns an array of size 1 x T, where T is the number of frames tracked
% in video expnum. 
%
% Access data for a fly indexed across all flies in all videos:
% flyi = 15
% expnum = data.summary.flies.expnum(flyi);
% flynum = data.summary.flies.flynum(flyi);
% y = data.exp(expnum).fly(flynum).y_mm;
% Returns an array of size 1 x T, where T is the number of frames tracked
% in video expnum. 
% 
% Access frames for all flies after the onset of the second period of
% activation:
% actidx = find(data.summary.activation.intervalnum == 2);
% deltat = 100;
% theta = zeros(0,deltat);
% for i = 1:numel(actidx),
%   inti = actidx(i);
%   expnum = data.summary.activation.expnum(inti);
%   t0 = data.summary.activation.startframe(inti);
%   for flynum = 1:numel(data.exp(expnum).fly),
%     theta = [theta;data.exp(expnum).fly(flynum).theta_rad(t0:t0+deltat-1)];
%   end
% end

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
      flycurr.(outnames_perframe{inti}) = nan(1,expinfo.nframes);
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
    