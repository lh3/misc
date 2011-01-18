// Compute genetic distance using the deCODE genetic map
if (arguments.length < 2) {
	print("Usage: k8 genetic_dist.js <decode.gmap> <reg.bed>");
	exit(1);
}

// read the genetic map
var f = new File(arguments[0]);
var $_, gmap = {};
var last = [null, '', 0, 0];
while (($_ = f.next()) != null) {
	var t = $_.split(/\s+/);
	if (t[3] == "cM") continue;
	if (last[0] != t[0]) {
		if (last[0]) gmap[last[0]].push([last[2], 1000000000, 0, 1]);
		last[2] = 0; gmap[t[0]] = [];
	}
	if (t[3] != "NA") gmap[t[0]].push([parseInt(last[2]), parseInt(t[2]), parseFloat(t[3])]);
	else gmap[t[0]].push([parseInt(last[2]), parseInt(t[2]), 0, 1]);
	last = t;
}
f.close();

// locate interval
function find_intv(a, x)
{
	var mid, min = 0, max = a.length;
	for (;;) { // a simplified binary search
		mid = Math.floor((min + max) / 2);
		if (x >= a[mid][0] && x <= a[mid][1]) break;
		if (x > a[mid][1]) min = mid + 1;
		else if (x < a[mid][0]) max = mid - 1;
	}
	return mid;
}

// read the BED
f = new File(arguments[1]);
while (($_ = f.next()) != null) {
	var t = $_.split(/\s+/);
	if (!gmap[t[0]]) continue; // chr not found
	var aa = gmap[t[0]];
	var beg = find_intv(aa, t[1]);
	if (t[2] > aa[aa.length-1][1]) t[2] = aa[aa.length-1][1];
	var end = find_intv(aa, t[2]);
	var gd = 0, pd = 0;
	for (i = beg + 1; i <= end - 1; ++i) {
		if (aa[i].length == 3) pd += aa[i][1] - aa[i][0];
		gd += aa[i][2];
	}
	if (beg == end) {
		if (aa[beg].length == 3) pd += end - beg;
		gd += (end - beg) / (aa[beg][1] - aa[beg][0]) * aa[beg][2]; 
	} else {
		if (aa[beg].length == 3) pd += aa[beg][1] - t[1];
		if (aa[end].length == 3) pd += t[2] - aa[end][0];
		gd += (aa[beg][1] - t[1]) / (aa[beg][1] - aa[beg][0]) * aa[beg][2];
		gd += (t[2] - aa[end][0]) / (aa[end][1] - aa[end][0]) * aa[end][2];
	}
	print(t[0]+"\t"+t[1]+"\t"+t[2]+"\t"+gd+"\t"+pd);
}
f.close();
