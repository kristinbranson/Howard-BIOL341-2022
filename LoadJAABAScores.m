% data = LoadJAABAScores(data,expdirs,scorefilestr,...)
%
% Load tracking information from experiment directories expdirs. 
% 
% Input:
% data: struct output of LoadTracking
% expdirs: cell of full paths to experiment directories. 
% scorefilestr: name of scores file
% 
% Output: 
% data: data is a struct. 
% We add data.exp(expnum).fly(flynum).(behaviorname) indicating the
% predictions at every frame. 

function data = LoadJAABAScores(data,expdirs,scorefilestr)

if ischar(expdirs),
  expdirs = {expdirs};
end

nexps = numel(expdirs);
for expnum = 1:nexps,
  sd = load(fullfile(expdirs{expnum},scorefilestr));
  behavior = sd.behaviorName;
  if iscell(behavior),
    behavior = behavior{1};
  end
  pred = sd.allScores.postprocessed;
  for flynum = 1:numel(data.exp(expnum).fly),
    data.exp(expnum).fly(flynum).(behavior) = pred{flynum};
  end
end

    