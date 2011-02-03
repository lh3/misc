matrix = {}

function matrix.T(a)
	local m, n, x = #a, #a[1], {};
	for i = 1, n do
		x[i] = {};
		for j = 1, m do x[i][j] = a[j][i] end
	end
	return x;
end

function matrix.mul(a, b)
	assert(#a[1] == #b);
	local m, n, p, x = #a, #a[1], #b[1], {};
	local c = matrix.T(b); -- transpose for efficiency
	for i = 1, m do
		x[i] = {}
		local xi = x[i];
		for j = 1, p do
			local sum, ai, cj = 0, a[i], c[j];
			for k = 1, n do sum = sum + ai[k] * cj[k] end
			xi[j] = sum;
		end
	end
	return x;
end
local n = arg[1] or 100;

function matgen()
	local a = {}
	for i = 1, n do
		a[i] = {}
		for j = 1, n do
			a[i][j] = math.random(-1, 1)
		end
	end
	return a;
end

local a = matrix.mul(matgen(), matgen());
