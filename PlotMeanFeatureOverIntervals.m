function hax = PlotMeanFeatureOverIntervals(data,featname,varargin)
% hax = PlotMeanFeatureOverIntervals(data,featname,...)
% Plots the mean feature featname around the start of each lights-on/off period.
% Inputs:
% data: struct with the data
% featname: field name for feature
% Output:
% hax: axes handles
% optional arguments:
% featlabel: string with label for y-axis (default: 'Feature (units)')
% minfeatplot: Lower limit for y-axis (default: 0)
% maxfeatplot: Upper limit for y-axis (default: []). If empty, will use the
% max of all data. 
% plotstderr: whether to plot the standard error (default: true)
% plotallflies: whether to plot individual flies (default: false)
% For the following names, if computing offset periods, replace 'on' with
% 'off'.
% onset: whether to plot the onset of activation (true) or the offset
% (false) (default: true). 
% meanfeatname: field name for mean of feature (default:
% ['mean_on_',featname])
% stderrfeatname: field name for stderr of feature (default:
% ['stderr_on_',featname])
% deltatname: field name to store interval information in (default:
% ['deltat_on_',featname]

nexps = numel(data.exp);

[featlabel,minfeatplot,maxfeatplot,plotstderr,...
  tfonset,meanfeatname,stderrfeatname,deltatname] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',0,'maxfeatplot',[],...
  'plotstderr',true,...
  'onset',true,...
  'meanfeatname','',...
  'stderrfeatname','',...
  'deltatname','');


if tfonset,
  onsetname = 'onset';
else
  onsetname = 'offset';
end

if isempty(meanfeatname)
  meanfeatname = sprintf('mean_%s_%s',onsetname,featname);
end
if isempty(stderrfeatname)
  stderrfeatname = sprintf('stderr_%s_%s',onsetname,featname);
end
if isempty(deltatname),
  deltatname = sprintf('deltat_%s_%s',onsetname,featname);
end

[exptypes,~,exptypeidx] = unique(data.summary.exps.type);

exptypecolors = lines(numel(exptypes));
ylim = [minfeatplot,maxfeatplot];
deltat = data.exp(1).stat.(deltatname);

naxc = 2;
naxr = ceil(nexps/naxc);

hax = gobjects(naxc,naxr);
for i = 1:nexps,

  fps = data.exp(i).summary.fps;
  expcolor = exptypecolors(exptypeidx(i),:);
  hax(i) = subplot(naxr,naxc,i);
  plot([0,0],ylim,'k-');
  hold on;

  meancurr = data.exp(i).stat.(meanfeatname);
 
  if plotstderr,
    stderrcurr = data.exp(i).stat.(stderrfeatname);

    plot((deltat(1):deltat(2))/fps,meancurr-stderrcurr,'-','Color',expcolor);
    plot((deltat(1):deltat(2))/fps,meancurr+stderrcurr,'-','Color',expcolor);
  end
  plot((deltat(1):deltat(2))/fps,meancurr,'-','Color',expcolor,'LineWidth',2);
  title(sprintf('%d: %s',i,data.exp(i).summary.type),'Interpreter','none');
end
ylabel(hax(1,naxr),featlabel);
xlabel(hax(1,naxr),'Time (s)');
fps = median(data.summary.exps.fps);
set(hax(1:nexps),'YLim',ylim,'XLim',deltat/fps);
linkaxes(hax(1:nexps));