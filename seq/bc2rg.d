#!/usr/bin/evn rdmd

import std.stdio, std.string, std.regex, klib;

void main(string[] args)
{
	if (args.length < 4) {
		writeln("Usage: rdmd bc2rg.d <sample> <barcode-list.txt> <aln.sam>");
		return;
	}
	auto fp = new ZFile(args[2]);
	ubyte[] l;
	string[string] bc2rg;
	char[] base = ['A', 'C', 'G', 'T', 'N'];
	string[] rg;
	while (fp.readto(l) >= 0) {
		auto t = std.string.split(cast(string)l);
		if (t.length == 0) continue;
		string r = args[1] ~ '-' ~ t[0] ~ '-' ~ t[1];
		rg ~= r;
		for (auto i = 0; i < t[0].length; ++i) {
			foreach (c; base) {
				char[] s;
				s.length = t[0].length;
				s[] = cast(char[])t[0];
				s[i] = c;
				bc2rg[cast(string)s] = r;
			}
		}
	}
	fp.close();
	fp = new ZFile(args[3]);
	bool first = true;
	auto re = regex(r"^\S+#([ACGTN]+)");
	while (fp.readto(l) >= 0) {
		if (l[0] != '@') {
			if (first) {
				foreach (p; rg)
					writefln("@RG\tID:%s\tSM:%s\tLB:%s", p, args[1], p);
				first = false;
			}
			auto c = match(cast(string)l, re).captures;
			if (c[1] in bc2rg) writeln(cast(char[])l, "\tRG:Z:", bc2rg[c[1]]);
			else writeln(cast(char[])l);
		} else writeln(cast(char[])l);
	}
}
