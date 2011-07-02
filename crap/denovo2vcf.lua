#!/usr/bin/env luajit

require('klib')

-- convert Conrad et al. (2011) validation table to VCF

for l in io.lines() do
	local t = l:split('\t')
	if t[1] ~= 'id' then
		local c = {t[15], t[16], t[11], t[12], t[13], t[14], t[39], t[40], 0, 0}
		for i = 17, 38, 2 do
			c[9], c[10] = c[9] + t[i], c[10] + t[i+1]
		end
		local tmp = t[5]:match('"(%S+)"')
		if tmp ~= nil then t[5] = tmp end
		if t[2] == '23' then t[2] = 'X' end
		local o = {t[2], t[3], t[1], t[4], t[5]}
		print(table.concat(o, '\t'), 30, '.',
			'Type='..t[50]..';FLT='..t[51]..';Set='..t[6]..t[7]..t[8]..';CpG='..t[49],
			'CUM:CSC',
			c[1]..','..c[2]..':'..t[41]..','..t[42],
			c[3]..','..c[4]..':'..t[43]..','..t[44],
			c[5]..','..c[6]..':'..t[45]..','..t[46],
			c[9]..','..c[10]..':'..t[47]..','..t[48],
			c[7]..','..c[8]..':0,0')
	else
		print('##fileformat=VCFv4.0')
		print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891\tNA12892\tChildren\tNA12877')
	end
end
