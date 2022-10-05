function m = parseExpName(expname)

m = regexp(expname,'^(?<exptype>.+)_(?<condition>.+)_Rig(?<rig>.+)_(?<timestamp>\d{8}T\d{6})$','names','once');
assert(~isempty(m));
