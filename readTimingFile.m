function [camerastart_datetime,protocolstart_datetime,...
  pulsestart_datetime,pulsestop_datetime] = readTimingFile(timinginfofile)

fid = fopen(timinginfofile,'r');
camerastart_datetime = [];
protocolstart_datetime = [];
pulsestart_datetime = [];
pulsestop_datetime = [];
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  if isempty(s),
    continue;
  end
  ss = regexp(s,',','split');
  for i = 1:2:numel(ss)-1,
    a = str2double(ss{1});
    b = str2double(ss{2});
    if isnan(b) && ~isnan(a),
      if strcmp(ss{2},'start_camera'),
        camerastart_datetime = a;
      elseif strcmp(ss{2},'run_experiment'),
        protocolstart_datetime = a;
      elseif strcmp(ss{2},'start_pulse'),
        pulsestart_datetime(end+1) = a; %#ok<AGROW> 
      elseif strcmp(ss{2},'stop_pulse'),
        pulsestop_datetime(end+1) = a; %#ok<AGROW> 
      end
    end
  end
end
fclose(fid);