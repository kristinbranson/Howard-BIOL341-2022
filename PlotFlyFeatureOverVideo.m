function hax = PlotFlyFeatureOverVideo(feat,activation)

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

nflies = numel(feat);
[featlabel,minfeatplot,maxfeatplot] = ...
  myparse(varargin,'featlabel','Feature (units)','minfeatplot',0,'maxfeatplot',[]);

