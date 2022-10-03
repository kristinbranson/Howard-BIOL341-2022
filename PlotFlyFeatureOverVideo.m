function hax = PlotFlyFeatureOverVideo(feat,fps,activation_startframes,activation_endframes,varargin)

% hax = PlotFlyFeatureOverVideo(feat,activation,fps,...)
% Inputs:
% feat: cell with an entry for each fly. feat{i} is the data feature for
% fly i and is a vector of size 1 x T
% activation: struct with information about activation periods 
% fps: frame rate of the camera
% Output:
% hax: axes handles
% optional arguments:
% featlabel: string with label for y-axis (default: 'Feature (units)')
% minfeatplot: Lower limit for y-axis (default: 0)
% maxfeatplot: Upper limit for y-axis (default: []). If empty, will use the
% max of all data. 

[nflies,T] = size(feat);

[featlabel,minfeatplot,maxfeatplot] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',0,'maxfeatplot',[]);
if isempty(maxfeatplot),
  maxfeatplot = max(feat(:));
end

colors = jet(nflies)*.7;
hax = gca;

isactivation = nargin >= 3 && ~isempty(activation_startframes);
hold(hax,'on');

for fly = 1:nflies,
  plot(hax,[0,(T+1)/fps],[fly-1,fly-1],'k:');
  if isactivation,
    acti = min(fly,size(activation_startframes,1));
    PlotActivationTimePatch(activation_startframes(acti,:),activation_endframes(acti,:),fps,[fly-1,fly],hax);
  end
  plot(hax,(1:T)/fps,fly-1+(feat(fly,:)-minfeatplot)/(maxfeatplot-minfeatplot),'-','Color',colors(fly,:));
end
 plot(hax,[0,(T+1)/fps],[nflies,nflies],'k:');
xlabel(hax,'Time (s)');
ylabel(hax,['Fly, ',featlabel]);
set(hax,'XLim',[0,(T+1)/fps],'YLim',[-.1,nflies+.1]);