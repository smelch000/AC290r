#!/bin/bash

[ $# -lt 1 ] || [ $# -gt 2 ] && echo "Usage: ${0##*/} <bgkflag file> <nodeownr file>" &&  exit 1

MYDIR=$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)

echo -n "Running change_bgk_layout.sh $1 $2..."
if [ $# -eq 1 ]; then
	$MYDIR/change_bgk_layout.sh "$1"
else
	$MYDIR/change_bgk_layout.sh "$1" "$2"
fi
echo "done"

echo "Running muphy..."
$MYDIR/parmuphy

exit $?
