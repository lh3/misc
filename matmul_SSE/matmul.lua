require("klua")

local n = arg[1] or 100;

-- generate matrix

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

local a = matgen();
local b = matgen();

local c = matrix.mul(a, b);
