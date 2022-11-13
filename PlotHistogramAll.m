function hax = PlotHistogramAll(data,featname,varargin)
% hax = PlotHistogramAll(data,featname,...)
% Plot the on and off histograms for feature featname

nexps = numel(data.exp);

[featlabel,fracfeatname,binsfeatname] = ...
  myparse(varargin,'featlabel','Feature (units)',...
  'fracfeatname','',...
  'binsfeatname','');

if isempty(fracfeatname)
  fracfeatname = sprintf('frac_all_%s',featname);
end
if isempty(binsfeatname)
  binsfeatname = sprintf('bincenters_all_%s',featname);
end

[exptypes,~,exptypeidx] = unique(data.summary.exps.type);

exptypecolors = lines(numel(exptypes));

naxc = min(nexps,2);
naxr = ceil(nexps/naxc);

hax = gobjects(naxc,naxr);
for i = 1:nexps,

  fps = data.exp(i).summary.fps;
  color = exptypecolors(exptypeidx(i),:);
  hax(i) = subplot(naxr,naxc,i);
  hold on;

  h = plot(hax(i),data.exp(i).hist.(binsfeatname),data.exp(i).hist.(fracfeatname),'-','Color',color,'LineWidth',2);
  set(hax(i),'XLim',data.exp(i).hist.(binsfeatname)([1,end]));

  title(sprintf('%d: %s',i,data.exp(i).summary.type),'Interpreter','none');
end
xlabel(hax(1,naxr),featlabel);
ylabel(hax(1,naxr),'Frac. frames');
linkaxes(hax);