function hax = PlotHistogramOnVsOff(data,featname,varargin)
% hax = PlotHistogramOnVsOff(data,featname,...)
% Plot the on and off histograms for feature featname

nexps = numel(data.exp);

[featlabel,onfracfeatname,onbinsfeatname,...
  offfracfeatname,offbinsfeatname] = ...
  myparse(varargin,'featlabel','Feature (units)',...
  'onfracfeatname','',...
  'onbinsfeatname','',...
  'offfracfeatname','',...
  'offbinsfeatname','');


if isempty(onfracfeatname)
  onfracfeatname = sprintf('frac_on_%s',featname);
end
if isempty(onbinsfeatname)
  onbinsfeatname = sprintf('bincenters_on_%s',featname);
end

if isempty(offfracfeatname)
  offfracfeatname = sprintf('frac_off_%s',featname);
end
if isempty(offbinsfeatname)
  offbinsfeatname = sprintf('bincenters_off_%s',featname);
end

[exptypes,~,exptypeidx] = unique(data.summary.exps.type);

exptypecolors = lines(numel(exptypes));
offcolor = [0,0,0];

naxc = min(nexps,2);
naxr = ceil(nexps/naxc);

hax = gobjects(naxc,naxr);
for i = 1:nexps,

  fps = data.exp(i).summary.fps;
  oncolor = exptypecolors(exptypeidx(i),:);
  hax(i) = subplot(naxr,naxc,i);
  hold on;

  assert(isequal(data.exp(i).hist.(onbinsfeatname),data.exp(i).hist.(offbinsfeatname)));
  hon = plot(hax(i),data.exp(i).hist.(onbinsfeatname),data.exp(i).hist.(onfracfeatname),'-','Color',oncolor,'LineWidth',2);
  hoff = plot(hax(i),data.exp(i).hist.(offbinsfeatname),data.exp(i).hist.(offfracfeatname),'-','Color',offcolor,'LineWidth',2);
  set(hax(i),'XLim',data.exp(i).hist.(onbinsfeatname)([1,end]));

  title(sprintf('%d: %s',i,data.exp(i).summary.type),'Interpreter','none');
  legend([hon,hoff],{'On','Off'});
end
xlabel(hax(1,naxr),featlabel);
ylabel(hax(1,naxr),'Frac. frames');
linkaxes(hax);