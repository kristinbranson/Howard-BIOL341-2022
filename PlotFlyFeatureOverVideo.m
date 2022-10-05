function hax = PlotFlyFeatureOverVideo(data,featname,varargin)

% hax = PlotFlyFeatureOverVideo(feat,activation,fps,...)
% Plot the features in matrix feat 
% Inputs:
% feat: matrix with a row for each fly. feat(i,:) is the data feature for
% fly i. nflies x T matrix.
% fps: frame rate of the camera
% activation_startframes: matrix with activation start periods for each
% fly. nflies x nactivation
% activation_endframes: matrix with activation end periods for each
% fly. nflies x nactivation
% Output:
% hax: axes handles
% Optional arguments:
% featlabel: string with label for y-axis (default: 'Feature (units)')
% minfeatplot: Lower limit for y-axis (default: 0)
% maxfeatplot: Upper limit for y-axis (default: []). If empty, will use the
% max of all data. 

nflies = size(data.summary.flies,1);
maxT = max(data.summary.exps.nframes./data.summary.exps.fps);

[featlabel,minfeatplot,maxfeatplot] = ...
  myparse(varargin,'featlabel',featname,'minfeatplot',[],'maxfeatplot',[]);
if isempty(maxfeatplot),
  maxfeatplot = minfeatplot;
  for flyi = 1:nflies,
    expnum = data.summary.flies.expnum(flyi);
    flynum = data.summary.flies.flynum(flyi);
    feat = data.exp(expnum).fly(flynum).(featname);
    maxfeatplot = max(maxfeatplot,max(feat(:)));
  end
end

colors = jet(nflies)*.7;
hax = gca;

hold(hax,'on');

for flyi = 1:nflies,
  expnum = data.summary.flies.expnum(flyi);
  flynum = data.summary.flies.flynum(flyi);
  plot(hax,[0,maxT],[flyi-1,flyi-1],'k:');
  fps = data.exp(expnum).summary.fps;
  T = data.exp(expnum).summary.nframes;
  if ~isempty(data.exp(expnum).activation),
    sf = data.exp(expnum).activation.startframe;
    ef = data.exp(expnum).activation.endframe;
    PlotActivationTimePatch(sf,ef,fps,[flyi-1,flyi],hax);
  end
  feat = data.exp(expnum).fly(flynum).(featname);
  plot(hax,(1:numel(feat))/fps,flyi-1+(feat-minfeatplot)/(maxfeatplot-minfeatplot),'-','Color',colors(flyi,:));
end
plot(hax,[0,maxT],[nflies,nflies],'k:');
xlabel(hax,'Time (s)','Interpreter','none');
ylabel(hax,['Fly, ',featlabel]);
set(hax,'XLim',[0,maxT],'YLim',[-.1,nflies+.1]);