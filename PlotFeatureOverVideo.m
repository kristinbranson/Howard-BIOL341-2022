function hax = PlotFeatureOverVideo(feat,meanfeat,activation,fps,varargin)

% hax = PlotFeatureOverVideo(feat,meanfeat,activation,fps,...)
% Inputs:
% feat: cell with an entry for each video. feat{i} is the data feature for
% video i and is a matrix of size T x ntraj. 
% meanfeat: for each video, mean over flies of the feature. meanfeat is a
% cell with an entry for each video, and meanfeat{i} is a T x 1 vector. 
% activation: struct with information about activation periods for each
% video
% fps: frame rate of the camera
% Output:
% hax: axes handles
% optional arguments:
% featlabel: string with label for y-axis (default: 'Feature (units)')
% minfeatplot: Lower limit for y-axis (default: 0)
% maxfeatplot: Upper limit for y-axis (default: []). If empty, will use the
% max of all data. 
% plotallflies: whether to plot individual flies (default: false)
% plotstderr: whether to plot the standard error (default: true)
% stderrfeat: for each video, standard error over flies of the feature.
% stderrfeat is a cell with an entry for each video, and meanfeat{i} is a T
% x 1 vector. This must be input if plotstderr is true. (default: {})
% genotypeidx: index indicating which genotype each video is from. array of
% size 1 x nvideos. (default: 1:nvideos)
% expnames: names of experiments. (default: {'','',...}).

nvideos = numel(feat);
[featlabel,minfeatplot,maxfeatplot,plotallflies,plotstderr,...
  stderrfeat,genotypeidx,expnames] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',0,'maxfeatplot',[],...
  'plotallflies',false,'plotstderr',true,...
  'stderrfeat',{},...
  'genotypeidx',1:nvideos,...
  'expnames',repmat({''},[1,nvideos]));

clf;
naxc = 2;
naxr = ceil(nvideos/naxc);
if isempty(maxfeatplot),
  maxfeatplot = max(cellfun(@(x) max(x(~isinf(x))),feat));
end
ylim = [minfeatplot,maxfeatplot];

genotypecolors = lines(max(genotypeidx));

% one set of axes per video
hax = gobjects(naxc,naxr);
for i = 1:nvideos,

  T = size(feat{i},1);  

  % create the axes
  hax(i) = subplot(naxr,naxc,i);
  hold(hax(i),'on');

  % plot when the lights are on
  PlotActivationTimePatch(activation.startframes{i},activation.endframes{i},fps,ylim,hax(i));

  % plot the speed for all the flies. this will create a different line for
  % each fly
  if plotallflies,
    plot(hax(i),(1:T)/fps,feat{i}','-');
    statcolor = 'k';
  else
    statcolor = genotypecolors(genotypeidx(i),:);
  end

  % plot standard error (standard deviation of mean)
  if plotstderr,
    stderrcurr = stderrfeat{i};
    plot(hax(i),(1:T)/fps,meanfeat{i}-stderrcurr,'-','Color',statcolor);
    plot(hax(i),(1:T)/fps,meanfeat{i}+stderrcurr,'-','Color',statcolor);
  end

  % plot the mean speed over all flies
  plot(hax(i),(1:T)/fps,meanfeat{i}','-','Linewidth',2,'Color',statcolor);

  title(hax(i),sprintf('Video %d %s',i,expnames{i}),'Interpreter','none');
  box(hax(i),'off');
  set(hax(i),'YLim',ylim,'XLim',[0,(T+1)/fps])
end
% all axes forced to have the same y limits
linkaxes(hax,'y');
xlabel(hax(1,end),'Time (s)');
ylabel(hax(1,end),featlabel);