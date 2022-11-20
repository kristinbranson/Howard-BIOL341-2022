function [fractime,countpos,countneg] = ComputeFractimeActivation(pred,t0s,t1s,t0_off,t1_off)

countpos = 0;
countneg = 0;
for i = 1:numel(t0s),
  if i == 1,
    tmin = 1;
  else
    tmin = t1s(i-1);
  end
  tmax = min(t1s(i),numel(pred));
  t0 = min(tmax,max(tmin,t0s(i)+t0_off));
  t1 = min(tmax,max(tmin,t1s(i)+t1_off));
  countpos = countpos + nnz(pred(t0:t1)>0);
  countneg = countneg + nnz(~isnan(pred(t0:t1)));
end
fractime = countpos / countneg;