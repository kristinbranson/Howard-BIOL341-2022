function [hax,noise] = PlotFracTime(data,fn,varargin)

[hax,noise] = myparse(varargin,'ax',[],'noise',[]);
if isempty(hax),
  hax = gca;
end
hold(hax,'on');

nexps = numel(data.exp);
nmale = nnz(data.summary.flies.sex=='m');
nfemale = nnz(data.summary.flies.sex=='f');
isboth = nmale > 0 && nfemale > 0;

isnoisein = ~isempty(noise);
if ~isnoisein,
  noise = cell(1,nexps);
end

[exptypes,~,exptypeidx] = unique(data.summary.exps.type);
exptypecolors = lines(numel(exptypes));

xlabels = [""];

xoff = 1;
for expnum = 1:nexps,

  exptype = data.summary.exps.type(expnum);
  color = exptypecolors(exptypeidx(expnum),:);

  fractime_exp = data.exp(expnum).stat.(fn);
  h = bar(hax,xoff,fractime_exp);
  set(h,'FaceColor',color*.5+.5);

  nflies = numel(data.exp(expnum).fly);
  if isnoisein,
    r = noise{expnum};
  else
    r = rand(1,nflies)*.2-.1;
    noise{expnum} = r;
  end
  fractime_fly = [data.exp(expnum).fly.(fn)];
  plot(hax,xoff+r,fractime_fly,'wo','MarkerFaceColor',color);
  xlabels(xoff) = exptype;

  xoff = xoff + 1;

  if isboth,
    fidx = [data.exp(expnum).fly.sex] == 'f';
    if any(fidx) && ~all(fidx),
      fn1 = [fn,'_female'];
      fractime_exp = data.exp(expnum).stat.(fn1);
      h = bar(hax,xoff,fractime_exp);
      set(h,'FaceColor',color*.75+.25);
      fidx = [data.exp(expnum).fly.sex] == 'f';
      plot(hax,xoff+r(fidx),fractime_fly(fidx),'+','Color',color,'LineWidth',2);
      xlabels(xoff) = exptype+", F";
      xoff = xoff + 1;


      fn1 = [fn,'_male'];
      fractime_exp = data.exp(expnum).stat.(fn1);
      h = bar(hax,xoff,fractime_exp);
      set(h,'FaceColor',color*.25+.75);
      plot(hax,xoff+r(~fidx),fractime_fly(~fidx),'x','Color',color,'LineWidth',2);
      xlabels(xoff) = exptype+", M";
      xoff = xoff + 1;
    end
  end
end
set(hax,'XTick',1:numel(xlabels),'XTickLabel',xlabels);
set(hax.XAxis,'TickLabelInterpreter','none');
set(hax,'XLim',[.5,xoff-.5]);
ylabel(hax,'Fraction of time')