require('klib')
require('bio')
local bit = require('bit')

local fp, last = io.xopen(arg[1]), ''
for l in fp:lines() do
	if l:sub(1, 1):find('%d') ~= nil then
		local t = l:split('\t', 12)
		t[4], t[5] = tonumber(t[4]), tonumber(t[5])
		if t[3] ~= last then
			print('>' .. t[3])
			last = t[3]
			collectgarbage('collect')
		end
		if t[6] == 'no-call' or t[6] == 'half' then
			print(string.rep('x', t[5] - t[4]))
		elseif t[6] == 'hom' or t[6] == 'hap' then
			if t[7] == 'ref' then print(string.rep('X', t[5] - t[4]))
			elseif (t[7] == 'sub' or t[7] == 'snp') and #t[9] == t[5] - t[4] then print(t[9])
			else print(string.rep('x', t[5] - t[4]))
			end
		elseif t[6] == 'het-ref' or t[6] == 'het-alt' then
			if (t[7] == 'sub' or t[7] == 'snp') and #t[9] == #t[10] and #t[9] == t[5] - t[4] then
				local s = {}
				for i = 1, #t[9] do
					local x = bit.bor(bio.nt16[t[9]:byte(i)], bio.nt16[t[10]:byte(i)])
					s[#s+1] = bio.ntrev:sub(x+1, x+1)
				end
				--print(t[4], t[9], t[10], table.concat(s))
				print(table.concat(s))
			else print(string.rep('x', t[5] - t[4]))
			end
		else print(string.rep('x', t[5] - t[4]))
		end
	end
end
