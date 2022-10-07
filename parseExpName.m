function m = parseExpName(expname,underscorein)

if nargin < 2,
  underscorein = 'condition';
end
if strcmpi(underscorein,'exptype'),
  pat = '^(?<exptype>[^_]+)_(?<condition>.+)';
elseif strcmpi(underscorein,'condition'),
  pat = '^(?<exptype>.+)_(?<condition>[^_]+)';
else
  error('Unknown underscorein = %s',underscorein);
end
pat = [pat,'_Rig(?<rig>.+)_(?<timestamp>\d{8}T\d{6})$'];
m = regexp(expname,pat,'names','once');
assert(~isempty(m));
m.type = [m.exptype,'_',m.condition];
