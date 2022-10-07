function data = CropTracking(data,expi,t0in,t1in)

T0 = data.exp(expi).summary.nframes;

if isfield(data.exp(expi).summary,'startframe') && data.exp(expi).summary.startframe > 1,
  t0 = t0in - data.exp(expi).summary.startframe + 1;
  t1 = t1in - data.exp(expi).summary.startframe + 1;
else
  t0 = t0in;
  t1 = t1in;
end
data.exp(expi).summary.startframe = t0in;
data.exp(expi).summary.endframe = t1in;
data.exp(expi).summary.nframes = t1in-t0in+1;
data.exp(expi).timestamps = data.exp(expi).timestamps(t0:t1);

data.summary.exps.startframe(expi) = data.exp(expi).summary.startframe;
data.summary.exps.endframe(expi) = data.exp(expi).summary.endframe;
data.summary.exps.nframes(expi) = data.exp(expi).summary.nframes;

% flies
nflies = numel(data.exp(expi).fly);
doremove = false(nflies,1);

fnsspecial = {'sex','startframe','endframe','nframes'};
skippedfns = {};
for flynum = 1:nflies,
  flycurr = data.exp(expi).fly(flynum);
  if flycurr.startframe > t1 || flycurr.endframe < t0,
    doremove(flynum) = true;
    continue;
  end
  flycurr.startframe = max(1,flycurr.startframe-t0+1);
  flycurr.endframe = min(t1,flycurr.endframe-t0+1);
  flycurr.nframes = flycurr.endframe-flycurr.startframe+1;
  flyi = find(data.summary.flies.expnum==expi & data.summary.flies.flynum==flynum);
  assert(numel(flyi)==1);
  data.summary.flies.startframe(flyi) = flycurr.startframe;
  data.summary.flies.endframe(flyi) = flycurr.endframe;
  data.summary.flies.nframes(flyi) = flycurr.nframes;

  fns = setdiff(fieldnames(flycurr),fnsspecial);

  for j = 1:numel(fns),
    fn = fns{j};
    off = T0-numel(flycurr.(fn));
    if abs(off) <= 3,
      flycurr.(fn) = flycurr.(fn)(t0:t1-off);
    else
      skippedfns = union(skippedfns,{fn});
    end
  end
  data.exp(expi).fly(flynum) = flycurr;
end
if any(doremove),
  data.exp(expi).fly(doremove) = [];
  flynumremove = find(doremove);
  flyremove = data.summary.flies.expnum==expi & ismember(data.summary.flies.flynum,flynumremove);
  data.summary.flies(flyremove,:) = [];
  data.summary.flies.flynum(data.summary.flies.expnum==expi) = 1:numel(data.exp(expi).fly);
  data.exp(expi).summary.nflies = numel(data.exp(expi).fly);
  data.summary.exps.nflies(expi) = data.exp(expi).summary.nflies;
end

% activation intervals
act = data.exp(expi).activation;
doremove = act.startframe > t1 | act.endframe < t0;
act(doremove,:) = [];
actremove = data.summary.activation.expnum == expi & ismember(data.summary.activation.intervalnum,find(doremove));
assert(nnz(actremove)==nnz(doremove));
data.summary.activation(actremove,:) = [];

if any(act.startframe < t0 | act.endframe > t1),
  warning('Cropping activation periods in the middle');
end
act.startframe = max(1,act.startframe-t0+1);
act.endframe = min(t1,act.endframe-t0+1);
data.exp(expi).activation = act;

oldintervalnum = find(~doremove);
for i = 1:size(act,1)
  acti = find(data.summary.activation.expnum==expi & data.summary.activation.intervalnum==oldintervalnum(i));
  assert(numel(acti)==1);
  data.summary.activation.intervalnum(acti) = i;
  data.summary.activation.startframe(acti) = act.startframe(i);
  data.summary.activation.endframe(acti) = act.endframe(i);
end
