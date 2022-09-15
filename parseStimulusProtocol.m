function [starttimes,endtimes,intensity,pulsewidth,pulseperiod,stepnums,stepstarts,isRGB] = ...
  parseStimulusProtocol(protocol,color)
% returns seconds

if nargin < 1,
  color = 'R';
end

isRGB = isfield(protocol,[color,'intensity']);

if isRGB,
  protocol = selectColor(protocol,color);
else
  protocol = modernize(protocol);
end

stepnums = [];
starttimes = [];
endtimes = [];
intensity = [];
pulsewidth = [];
pulseperiod = [];
stepstarts = nan(1,numel(protocol.stepNum));
stepstart = 0; % ms
for stepi = 1:numel(protocol.stepNum), 
  stepstarts(stepi) = stepstart;
  niter = protocol.iteration(stepi);
  ledonstart = protocol.delayTime(stepi)*1000 + stepstart; % delay time in s
  ondur = protocol.pulsePeriod(stepi)*protocol.pulseNum(stepi);
  offdur = protocol.offTime(stepi);
  totaldur = ondur+offdur;
  starttimes_step = ledonstart+totaldur*(0:niter-1);
  endtimes_step = starttimes_step+ondur;
  stepstart = stepstart + protocol.duration(stepi);
  starttimes = [starttimes,starttimes_step]; %#ok<AGROW> 
  endtimes = [endtimes,endtimes_step]; %#ok<AGROW> 
  intensity = [intensity,protocol.intensity(stepi)+zeros(1,niter)]; %#ok<AGROW> 
  pulseperiod = [pulseperiod,protocol.pulsePeriod(stepi)+zeros(1,niter)]; %#ok<AGROW> 
  pulsewidth = [pulsewidth,protocol.pulseWidth(stepi)+zeros(1,niter)]; %#ok<AGROW> 
  stepnums = [stepnums,stepi+zeros(1,niter)]; %#ok<AGROW> 
end

% convert to seconds
starttimes = starttimes / 1000;
endtimes = endtimes / 1000;
stepstarts = stepstarts / 1000;

function protocolout = selectColor(protocolin,color)

protocolout = struct;

fnscopy = {'stepNum','duration','delayTime','ProtocolData','ProtocolHeader'};
fnscolor = {'intensity','pulseWidth','pulsePeriod','pulseNum','offTime','iteration'};
for i = 1:numel(fnscopy),
  fn = fnscopy{i};
  protocolout.(fn) = protocolin.(fn);
end
for i = 1:numel(fnscolor),
  fn = fnscolor{i};
  protocolout.(fn) = protocolin.([color,fn]);
end

function protocolout = modernize(protocolin)

protocolout = struct;
fnscopy = {'stepNum','duration','delayTime','ProtocolData','ProtocolHeader','intensity','pulseNum','offTime','iteration'};
for i = 1:numel(fnscopy),
  fn = fnscopy{i};
  protocolout.(fn) = protocolin.(fn);
end
protocolout.pulseWidth = protocolin.pulseWidthSP;
protocolout.pulsePeriod = protocolin.pulsePeriodSP;
