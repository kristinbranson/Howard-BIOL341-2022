function hax = PlotFeatureOverTime(feat,activation,fps,varargin)

% optional arguments:
% featlabel: string with label for y-axis (default: 'Feature (units)')
% minfeatplot: Lower limit for y-axis (default: 0)
% maxfeatplot: Upper limit for y-axis (default: []). If empty, will use the
% max of all data. 
% plotallflies: whether to plot individual flies (default: false)
% plotstderr: whether to plot the standard error (default: true)

nvideos = numel(feat);
[featlabel,minfeatplot,maxfeatplot,plotallflies,plotstderr,...
  meanfeatpervideo,stderrfeatpervideo,genotypeidx,expnames] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',0,'maxfeatplot',[],...
  'plotallflies',false,'plotstderr',true,...
  'meanfeat',{},'stderrfeat',{},...
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
hax = gobjects(nvideos,1);
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
    stderrcurr = stderrfeatpervideo{i};
    plot(hax(i),(1:T)/fps,meanfeatpervideo{i}-stderrcurr,'-','Color',statcolor);
    plot(hax(i),(1:T)/fps,meanfeatpervideo{i}+stderrcurr,'-','Color',statcolor);
  end

  % plot the mean speed over all flies
  plot(hax(i),(1:T)/fps,meanfeatpervideo{i}','-','Linewidth',2,'Color',statcolor);

  title(hax(i),sprintf('Video %d %s',i,expnames{i}),'Interpreter','none');
  box(hax(i),'off');
  set(hax(i),'YLim',ylim,'XLim',[0,(T+1)/fps])
end
% all axes forced to have the same y limits
linkaxes(hax,'y');
xlabel(hax(end),'Time (s)');
ylabel(hax(end),featlabel);