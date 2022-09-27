function PlotFeatureOverIntervals(feat,meanfeat,activation,fps,varargin)

% hax = PlotFeatureOverIntervals(feat,meanfeat,activation,fps,...)
% Inputs:
% feat: cell with an entry for each video. feat{i} is the data feature for
% video i and is a matrix of size T x ntraj. 
% meanfeat: for each video, mean over flies of the feature. meanfeat is a
% cell with an entry for each video, and meanfeat{i} is a T x 1 vector. 
% activation: struct with information about activation periods for each
% video
% fps: frame rate of the camera
% Output:
% hax: axes handles
% optional arguments:
% featlabel: string with label for y-axis (default: 'Feature (units)')
% minfeatplot: Lower limit for y-axis (default: 0)
% maxfeatplot: Upper limit for y-axis (default: []). If empty, will use the
% max of all data. 
% plotallflies: whether to plot individual flies (default: false)
% plotstderr: whether to plot the standard error (default: true)
% stderrfeat: for each video, standard error over flies of the feature.
% stderrfeat is a cell with an entry for each video, and meanfeat{i} is a T
% x 1 vector. This must be input if plotstderr is true. (default: {})
% genotypeidx: index indicating which genotype each video is from. array of
% size 1 x nvideos. (default: 1:nvideos)
% nframespre_plot: how much before the lights on to plot (default: 1*fps)
% nframespost_plot: how much after the lights on to plot (default: 5*fps)
% nstimuliplot: how many stimulus periods to plot. if there are less than
% nstimuliplot periods, all stimuli will be plotted (default: inf). 
% plotonset: whether to plot the onset of activation (true) or the offset
% (false) (default: true). 

nvideos = numel(feat);

[featlabel,minfeatplot,maxfeatplot,plotallflies,plotstderr,...
  stderrfeat,genotypeidx,nframespre_plot,nframespost_plot,...
  nstimuliplot,plotonset,expnames] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',0,'maxfeatplot',[],...
  'plotallflies',false,'plotstderr',true,...
  'stderrfeat',{},...
  'genotypeidx',1:nvideos,...
  'nframespre_plot',1*fps,...
  'nframespost_plot',5*fps,...
  'nstimuliplot',inf,...
  'plotonset',true,...
  'expnames',repmat({''},[1,nvideos]));

ylim = [minfeatplot,maxfeatplot];

if plotonset,
  allts = activation.startframes;
else
  allts = activation.endframes;
end
nstimuli = nan(1,nvideos);
for i = 1:nvideos,
  nstimuli(i) = numel(allts{i});
end
nstimuliplot = min(nstimuliplot,max(nstimuli));
middlecol = ceil(nstimuliplot/2);

genotypecolors = lines(max(genotypeidx));

hax = gobjects(nvideos,nstimuliplot);

for i = 1:nvideos,
  T = size(feat{i},1);
  for j = 1:nstimuliplot,
    if j > nstimuli(i),
      continue;
    end

    % which axes to plot in
    axi = sub2ind([nstimuliplot,nvideos],j,i);
    hax(i,j) = subplot(nvideos,nstimuliplot,axi);

    % which frames to plot
    t = allts{i}(j);

    t0 = max(1,t-nframespre_plot);
    t1 = min(T,t+nframespost_plot);

    hold(hax(i,j),'on');
    % plot each fly
    if plotallflies,
      plot(hax(i,j),(t0-t:t1-t)/fps,feat{i}(t0:t1,:)','-');
      statcolor = 'k';
    else
      statcolor = genotypecolors(genotypeidx(i),:);
    end
    if plotstderr,
      stderrcurr = stderrfeat{i}(t0:t1);
      plot(hax(i,j),(t0-t:t1-t)/fps,meanfeat{i}(t0:t1)-stderrcurr,'-','Color',statcolor);
      plot(hax(i,j),(t0-t:t1-t)/fps,meanfeat{i}(t0:t1)+stderrcurr,'-','Color',statcolor);
    end

    % plot a vertical line when the lights turn on
    plot(hax(i,j),[0,0],ylim,'r-');
    % plot the average speed of all flies
    if plotallflies,
      lw = 2;
    else
      lw = 1;
    end
    plot(hax(i,j),(t0-t:t1-t)/fps,meanfeat{i}(t0:t1)','-','LineWidth',lw,'color',statcolor);

    % only show which stimulus this is in the top row of axes
    if i == 1,
      title(hax(i,j),sprintf('Stimulus %d',j),'Interpreter','none');
    end

    if j == middlecol,
      xlabel(hax(i,j),expnames{i},'Interpreter','none');
    end
    box(hax(i,j),'off');
  end
  % only show which video this is in the left column of axes
  ylabel(hax(i,1),sprintf('Video %d',i));
end

% force axes to share the same x and y limits
linkaxes(hax(ishandle(hax)));
ylabel(hax(end,1),{featlabel,sprintf('Video %d',nvideos)});
xlabel(hax(end,1),'Time (s)');
set(hax(ishandle(hax)),'YLim',ylim,'XLim',[-nframespre_plot,nframespost_plot]/fps);