function hax = PlotMeanFeatureOverIntervals(meanfeat,stderrfeat,fps,nframespre,nframespost,varargin)

nvideos = size(meanfeat,2);

[featlabel,minfeatplot,maxfeatplot,plotstderr,...
  genotypeidx,expnames] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',0,'maxfeatplot',[],...
  'plotstderr',true,...
  'genotypeidx',1:nvideos,...
  'expnames',{});

if isempty(expnames),
  expnames = cell(1,nvideos);
  for i = 1:nvideos,
    expnames{i} = sprintf('Video %d',i);
  end
end

genotypecolors = lines(max(genotypeidx));
ylim = [minfeatplot,maxfeatplot];

naxc = 2;
naxr = ceil(nvideos/naxc);

hax = gobjects(naxc,naxr);
for i = 1:nvideos,

  videocolor = genotypecolors(genotypeidx(i),:);
  hax(i) = subplot(naxr,naxc,i);
  plot([0,0],ylim,'k-');
  hold on;

  if plotstderr,
    plot((-nframespre:nframespost)/fps,meanfeat(:,i)-stderrfeat(:,i),'-','Color',videocolor);
    plot((-nframespre:nframespost)/fps,meanfeat(:,i)+stderrfeat(:,i),'-','Color',videocolor);
  end
  plot((-nframespre:nframespost)/fps,meanfeat(:,i),'-','Color',videocolor,'LineWidth',2);
  title(expnames{i},'Interpreter','none');
end
ylabel(hax(1,naxr),featlabel);
xlabel(hax(1,naxr),'Time (s)');
set(hax,'YLim',ylim,'XLim',[-nframespre,nframespost]/fps);
linkaxes(hax);