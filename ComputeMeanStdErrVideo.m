function data = ComputeMeanStdErrVideo(data,featname,meanfeatname,stdfeatname,stderrfeatname)
% data = ComputeMeanStdErrVideo(data,featname,meanfeatname,stdfeatname,stderrfeatname)
% For each video, computes the mean, standard deviation, and standard error
% over flies for the feature featname. The outputs are stored in fields meanfeatname,
% stdfeatname, and stderrfeatname. If meanfeatname, stdfeatname, and
% stderrfeatname are not specified, they are set to ['mean_',featname],
% ['std_',featname], and ['stderr_',featname], respectively. 

if nargin < 3,
  meanfeatname = ['mean_',featname];
end
if nargin < 4,
  stdfeatname = ['std_',featname];
end
if nargin < 5,
  stderrfeatname = ['stderr_',featname];
end

nexps = numel(data.exp);

for expnum = 1:nexps,
  T = numel(data.exp(expnum).fly(1).(featname));

  % average over flies
  meanfeat = zeros(1,T);
  ncurr = zeros(1,T);
  for flynum = 1:numel(data.exp(expnum).fly),
    idx = ~isnan(data.exp(expnum).fly(flynum).(featname));
    meanfeat(idx) = meanfeat(idx) + data.exp(expnum).fly(flynum).(featname)(idx);
    ncurr(idx) = ncurr(idx) + 1;
  end
  meanfeat = meanfeat ./ ncurr;
  meanfeat(ncurr==0) = nan;
  data.exp(expnum).stat.(meanfeatname) = meanfeat;

  % standard deviation over flies
  stdfeat = zeros(1,T);
  for flynum = 1:numel(data.exp(expnum).fly),
    idx = ~isnan(data.exp(expnum).fly(flynum).(featname));
    stdfeat(idx) = stdfeat(idx) + (data.exp(expnum).fly(flynum).(featname)(idx)-meanfeat(idx)).^2;
  end
  stdfeat = sqrt(stdfeat ./ ncurr);
  stdfeat(ncurr==0) = nan;
  data.exp(expnum).stat.(stdfeatname) = stdfeat;

  % standard error over flies
  data.exp(expnum).stat.(stderrfeatname) = stdfeat./sqrt(ncurr);
  data.exp(expnum).stat.(stderrfeatname)(ncurr==0) = nan;
  
end