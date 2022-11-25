{ pkgs }: {
	deps = [
		pkgs.gnuplot
  pkgs.clang_12
		pkgs.ccls
		pkgs.gdb
		pkgs.gnumake
	];
}