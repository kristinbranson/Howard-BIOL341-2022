% [h,hwing] = DrawWingedFly(x,y,theta,a,b,xwingl,ywingl,xwingr,ywingr,...)
%
% Draw a fly's current pose and position.
%
% Inputs: 
% x: x-coordinate of centroid
% y: y-coordinate of centroid
% theta: orientation in radians
% a: quarter-major axis length
% b: quarter-minor axis length
% xwingl: x-coordinate of left wing tip
% ywingl: y-coordinate of left wing tip
% xwingr: x-coordinate of right wing tip
% ywingr: y-coordinate of right wing tip
%
% Outputs:
% h: handle to body line
% hwing: handle to wing line
% 
% Optional arguments:
% 
% style: 'ellipse' or 'triangle'. Whether to draw an ellipse or triangle.
%
% Any additional arguments are passed to plot commands.

function [h,hwing] = DrawWingedFly(x,y,theta,a,b,xwingl,ywingl,xwingr,ywingr,varargin)

hwing = [];
[style,leftovers] = myparse_nocheck(varargin,'style','triangle');

if strcmpi(style,'ellipse'),

  h = ellipsedraw(a*2,b*2,x,y,theta,leftovers{:});

else

  % isosceles triangle not yet rotated or centered
  pts = [-a*2,-b*2
    -a*2,b*2
    a*2,0];
  
  % rotate
  costheta = cos(theta);
  sintheta = sin(theta);
  R = [costheta,sintheta;-sintheta,costheta];
  pts = pts*R;

  % translate
  pts(:,1) = pts(:,1) + x;
  pts(:,2) = pts(:,2) + y;

  % plot
  h = plot(pts([1:3,1],1),pts([1:3,1],2),leftovers{:});

end

if ~exist('xwingl','var'),
  return;
end

hwing = plot([xwingl,x,xwingr],[ywingl,y,ywingr],'-',leftovers{:});