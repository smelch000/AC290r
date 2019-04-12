#!/bin/bash

AWK=awk
AWK=gawk # for mac

[ $# -lt 1 ] || [ $# -gt 2 ] && echo "Usage: ${0##*/} <bgkflag file> <nodeownr file>" &&  exit 1

[ ! -f "$1" ] && echo "Cannot read $1!" && exit 1
[ $# -eq 2 ] && [ ! -f "$2" ] && "Cannot read $2!" && exit 1

BGKFILE="$1"
if [ $# -eq 2 ]; then
	PASTE_CMD="paste - \"$2\""
else
	PASTE_CMD="sed -e 's/$/ 0 0/'"
fi

HEADER=$(head -n 1 "$BGKFILE" | tr -s ' ' | sed -e 's/^[ ]\+//')

tail -n +3 "$BGKFILE" | eval $PASTE_CMD |\
$AWK -v header="$HEADER" '
{
	if (!seen[$6]) {
		seen[$6]=1
		print header > "bgkflag_"$6".hdr"
		nodes[$6,1] = nodes[$6,2] = nodes[$6,3] = nodes[$6,4] = 0
		print $1" "$2" "$3" "$4 > "bgkflag_"$6".dat"
	} else
		print $1" "$2" "$3" "$4 >> "bgkflag_"$6".dat"
	nodes[$6,$4]++
}
END {
	for(proc in seen)
		print "5 "nodes[proc,1]" "nodes[proc,2]" "nodes[proc,3]" "nodes[proc,4] >> "bgkflag_"proc".hdr"
}'

if [ $# -eq 1 ]; then
	mv bgkflag_0.dat bgkflag.dat
	mv bgkflag_0.hdr bgkflag.hdr
fi
