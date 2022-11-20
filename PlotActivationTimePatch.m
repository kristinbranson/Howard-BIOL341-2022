function h = PlotActivationTimePatch(t0s,t1s,fps,ylim,hax)
% PlotActivationTimePatch(t0s,t1s,fps,ylim,hax)
% Plots a pink box over windows between t0s(i) and t1s(i) for each i, with
% box heights defined by ylim. 

lightcolor = [1,.7,.7]; % pink
startcolor = [1,0,0];
startlw = 2;

h = gobjects(numel(t0s),1);
for j = 1:numel(t0s),
  t0 = t0s(j);
  t1 = t1s(j);
  h(j) = patch(hax,[t0,t0,t1,t1,t0]/fps,ylim([1,2,2,1,1]),lightcolor,'LineStyle','none');
end

% plot starts of activation periods
hline = gobjects(numel(t0s),1);
for j = 1:numel(t0s),
  t0 = t0s(j);
  hline(j) = plot(hax,[t0,t0]/fps,ylim,'-','Color',startcolor,'LineWidth',startlw);
end