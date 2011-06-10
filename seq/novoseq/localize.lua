#!/usr/bin/env luajit

-- localize unmapped segments

require('klib')
bit = require('bit')

local min_mq = 30

for o, a in os.getopt(arg, 'q:') do
	if o == 'q' then min_mq = a end
end

if #arg == 0 then
	print("Usage: localize.lua [-q minmq] <in.list> <in.sam>")
	os.exit(1)
end

local MAX_LEN = 1000000000

-- read the list
local h = {} -- the hash table
local fp = io.xopen(arg[1])
for l in fp:lines() do
	local t = l:split('\t')
	local s, b, e = t[1]:match('^(%S+):(%d+)-(%d+)')
	b, e = tonumber(b) - 1, tonumber(e)
	if h[s] == nil then
		h[s] = {{b, e, -1, MAX_LEN, 0, 0, '*', -1, '*', -1}} -- left, right, anchor_left, anchor_right, lqual, rqual, lchr, lpos, rchr, rpos
	else
		table.insert(h[s], {b, e, -1, MAX_LEN, 0, 0, '*', -1, '*', -1})
	end
end
fp:close()

-- parse the SAM
fp = io.xopen(arg[2])
local last = ''
for l in fp:lines() do
	local t = l:split('\t', 8)
	local ht = h[t[1]]
	if ht ~= nil then
		t[2], t[4], t[5] = tonumber(t[2]), tonumber(t[4]), tonumber(t[5])
		if t[5] >= min_mq then
			local qb, qe, tb, te
			local x, y, c = t[4] - 1, 0, {0, 0}
			for j, o in t[6]:gmatch('(%d+)(%S)') do
				j = tonumber(j)
				if o == 'M' then x, y = x + j, y + j
				elseif o == 'I' then y = y + j
				elseif o == 'D' then x = x + j
				 elseif o == 'H' or o == 'S' then
					c[y == 0 and 1 or 2] = j
					y = y + j
				end
			end
			local strand = bit.band(t[2], 0x10) == 0x10 and -1 or 1
			if strand > 0 then -- forward strand
				qb, qe, tb, te = c[1], y - c[2], t[4] - 1, x
			else
				qb, qe, tb, te = c[2], y - c[1], x, t[4] - 1
			end
			--print(t[1], qb, qe, t[3], tb, te, (strand > 0 and '+' or '-'), t[5])
			for _, v in pairs(ht) do
				if qe <= v[1] then
					if qe > v[3] then
						v[3], v[5], v[7], v[8] = qe, t[5], t[3], te
					end
				elseif qb >= v[2] then
					if qb < v[4] then
						v[4], v[6], v[9], v[10] = qb, t[5], t[3], tb
					end
				end
			end
		end
	end
	if t[1] ~= last and last ~= '' and h[last] ~= nil then
		for _, v in pairs(h[last]) do
			if v[4] == MAX_LEN then v[4] = -1 end
			print(last, table.concat(v, '\t'))
		end
	end
	last = t[1]
end
-- output the last if there is one
if last ~= '' and h[last] ~= nil then
	for _, v in pairs(h[last]) do
		if v[4] == MAX_LEN then v[4] = -1 end
		print(last, table.concat(v, '\t'))
	end
end
