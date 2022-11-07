function hax = PlotFeatureOverVideo(data,featname,varargin)

% hax = PlotFeatureOverVideo(data,featname,...)
% Inputs:
% data: struct with the data
% featname: field name for feature
% Output:
% hax: axes handles
% Optional arguments:
% featlabel: string with label for y-axis (default: 'Feature (units)')
% minfeatplot: Lower limit for y-axis (default: 0)
% maxfeatplot: Upper limit for y-axis (default: []). If empty, will use the
% max of all data. 
% plotallflies: whether to plot individual flies (default: false)
% plotstderr: whether to plot the standard error (default: true)
% meanfeatname: field name for mean of feature (default:
% ['mean_',featname])
% stderrfeatname: field name for stderr of feature (default:
% ['stderr_',featname])

nexps = numel(data.exp);
[featlabel,minfeatplot,maxfeatplot,plotallflies,plotstderr,...
  meanfeatname,stderrfeatname] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',[],'maxfeatplot',[],...
  'plotallflies',false,'plotstderr',true,...
  'meanfeatname',['mean_',featname],...
  'stdferreatname',['stderr_',featname]);

clf;
naxc = 2;
naxr = ceil(nexps/naxc);
if isempty(maxfeatplot),
  maxfeatplot = FeaturePrctiles(data,featname,100);
end
if isempty(minfeatplot),
  minfeatplot = FeaturePrctiles(data,featname,0);
end

ylim = [minfeatplot,maxfeatplot];

[exptypes,~,exptypeidx] = unique(data.summary.exps.type);
exptypecolors = lines(numel(exptypes));

% one set of axes per video
hax = gobjects(naxc,naxr);
for i = 1:nexps,

  T = numel(data.exp(i).fly(1).(featname));
  fps = data.exp(i).summary.fps;

  % create the axes
  hax(i) = subplot(naxr,naxc,i);
  hold(hax(i),'on');

  % plot when the lights are on
  PlotActivationTimePatch(data.exp(i).activation.startframe,data.exp(i).activation.endframe,fps,ylim,hax(i));

  % plot the speed for all the flies. this will create a different line for
  % each fly
  
  if plotallflies,
    flycolors = jet(numel(data.exp(i).fly))*.7;
    for j = 1:numel(data.exp(i).fly),
      plot(hax(i),(1:T)/fps,data.exp(i).fly(j).feat','-','Color',flycolors(j,:));
    end
    statcolor = 'k';
  else
    statcolor = exptypecolors(exptypeidx(i),:);
  end

  meancurr = data.exp(i).stat.(meanfeatname);

  % plot standard error (standard deviation of mean)
  if plotstderr,
    stderrcurr = data.exp(i).stat.(stderrfeatname);
    plot(hax(i),(1:T)/fps,meancurr-stderrcurr,'-','Color',statcolor);
    plot(hax(i),(1:T)/fps,meancurr+stderrcurr,'-','Color',statcolor);
  end

  % plot the mean speed over all flies
  plot(hax(i),(1:T)/fps,meancurr,'-','Linewidth',2,'Color',statcolor);

  title(hax(i),sprintf('Video %d %s',i,data.exp(i).summary.type),'Interpreter','none');
  box(hax(i),'off');
  set(hax(i),'YLim',ylim,'XLim',[0,(T+1)/fps])
end
% all axes forced to have the same y limits
linkaxes(hax,'y');
xlabel(hax(1,end),'Time (s)');
ylabel(hax(1,end),featlabel);