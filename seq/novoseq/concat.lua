require('klib')
require('bio')

local fp = io.xopen(arg[1])

local nlen, name, seqs, bed, coor = 20, 'hs37d5', {}, {}, 0
for n, s in bio.readseq(fp) do
	if #seqs ~= 0 then coor, seqs[#seqs+1] = coor + nlen, string.rep('n', nlen) end
	seqs[#seqs+1] = s
	bed[#bed+1] = {coor, coor + #s, n}
	coor = coor + #s
end

print('FA\t>'..name)
local s = table.concat(seqs)
for l = 1, #s, 60 do
	print('FA\t'..s:sub(l, l + 60 - 1))
end
for _, v in ipairs(bed) do
	local t = {'BD', name, v[1], v[2], v[3]}
	print(table.concat(t, '\t'))
end
