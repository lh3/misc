#!/usr/bin/env luajit

require('klib')

if #arg == 0 then
	print("Usage: luajit cg2vcf.lua <masterVarBeta.tsv.bz2>")
	os.exit(1)
end

local fp = io.xopen(arg[1])

for l in fp:lines() do
	local t = l:split("\t", 14)
	if t[1]:match('%d') and t[6] ~= 'no-call' and t[7] ~= 'ref' then
		if #t[8] == #t[9] and t[10] == '?' then -- half substitution
			for i = 1, #t[8] do
				if t[8]:byte(i) ~= t[9]:byte(i) then
					local info = {'SubLen='..#t[8], 'HALF'}
					print(t[3], t[4] + i, '.', t[8]:sub(i, i), t[9]:sub(i, i), t[11], '.', table.concat(info, ';'))
				end
			end
		elseif #t[8] == #t[9] and #t[8] == #t[10] then -- full substitution
			for i = 1, #t[8] do
				if t[8]:byte(i) ~= t[9]:byte(i) or t[8]:byte(i) ~= t[10]:byte(i) then -- snp
					local ref, a1, a2, alt, qual = t[8]:sub(i, i), t[9]:sub(i, i), t[10]:sub(i, i)
					info = {'SubLen='..#t[8], a1 == a2 and 'HOM' or 'HET'}
					if ref ~= a1 then alt = a1 end
					if ref ~= a2 and a1 ~= a2 then alt = alt and alt .. ',' .. a2 or a2 end
					if ref ~= a1 and ref ~= a2 then qual = tonumber(t[11]) > tonumber(t[12]) and t[11] or t[12] -- take the larger
					else qual = ref ~= a1 and t[11] or t[12] end -- take the VAF of the variant allele
					print(t[3], t[4] + i, '.', ref, alt, qual, '.', table.concat(info, ';'))
				end
			end
		elseif t[10] ~= '?' then -- full indel
			local info, alt, qual = {'INDEL', t[9] == t[10] and 'HOM' or 'HET'}
			if t[8] ~= t[9] then alt = 'N'..t[9] end
			if t[8] ~= t[10] and t[9] ~= t[10] then alt = alt and alt .. ',N' .. t[10] or 'N'..t[10] end
			if t[8] ~= t[9] and t[8] ~= t[10] then qual = tonumber(t[11]) > tonumber(t[12]) and t[11] or t[12]
			else qual = t[8] ~= t[9] and t[11] or t[12] end
			print(t[3], t[4], '.', 'N'..t[8], alt, qual, '.', table.concat(info, ';'))
		elseif t[10] == '?' then -- half indel
			local info = {'INDEL', 'HALF'}
			print(t[3], t[4], '.', 'N'..t[8], 'N'..t[9], t[11], '.', table.concat(info, ';'))
		else
			print(l)
		end
	end
end

fp:close()
