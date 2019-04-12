#!/bin/bash

[ $# -ne 6 ] && echo "Usage: ${0##*/} <x center> <y center> <z center> <x size> <y size> <z size>" &&  exit 1

MINX=$(echo "scale=4; $1-($4/2.00)" | bc)
MAXX=$(echo "scale=4; $1+($4/2.00)" | bc)

MINY=$(echo "scale=4; $2-($5/2.00)" | bc)
MAXY=$(echo "scale=4; $2+($5/2.00)" | bc)

MINZ=$(echo "scale=4; $3-($6/2.00)" | bc)
MAXZ=$(echo "scale=4; $3+($6/2.00)" | bc)

cat <<-EOF
solid CUBE
facet normal -1.0 0.0 0.0
outer loop
vertex $MINX $MINY $MINZ
vertex $MINX $MINY $MAXZ
vertex $MINX $MAXY $MINZ
endloop
endfacet
facet normal -1.0 0.0 0.0
outer loop
vertex $MINX $MINY $MAXZ
vertex $MINX $MAXY $MINZ
vertex $MINX $MAXY $MAXZ
endloop
endfacet
facet normal 0.0 -1.0 0.0
outer loop
vertex $MINX $MINY $MAXZ
vertex $MAXX $MINY $MINZ
vertex $MAXX $MINY $MAXZ
endloop
endfacet
facet normal 0.0 -1.0 0.0
outer loop
vertex $MINX $MINY $MAXZ
vertex $MINX $MINY $MINZ
vertex $MAXX $MINY $MINZ
endloop
endfacet
facet normal 0.0 1.0 0.0
outer loop
vertex $MAXX $MAXY $MAXZ
vertex $MAXX $MAXY $MINZ
vertex $MINX $MAXY $MAXZ
endloop
endfacet
facet normal 0.0 1.0 0.0
outer loop
vertex $MINX $MAXY $MAXZ
vertex $MAXX $MAXY $MINZ
vertex $MINX $MAXY $MINZ
endloop
endfacet
facet normal 0.0 0.0 -1.0
outer loop
vertex $MINX $MINY $MINZ
vertex $MAXX $MINY $MINZ
vertex $MINX $MAXY $MINZ
endloop
endfacet
facet normal 0.0 0.0 -1.0
outer loop
vertex $MINX $MAXY $MINZ
vertex $MAXX $MINY $MINZ
vertex $MAXX $MAXY $MINZ
endloop
endfacet
facet normal 0.0 0.0 1.0
outer loop
vertex $MINX $MINY $MAXZ
vertex $MAXX $MINY $MAXZ
vertex $MINX $MAXY $MAXZ
endloop
endfacet
facet normal 0.0 0.0 1.0
outer loop
vertex $MINX $MAXY $MAXZ
vertex $MAXX $MINY $MAXZ
vertex $MAXX $MAXY $MAXZ
endloop
endfacet
facet normal 1.0 0.0 0.0
outer loop
vertex $MAXX $MINY $MAXZ
vertex $MAXX $MINY $MINZ
vertex $MAXX $MAXY $MINZ
endloop
endfacet
facet normal 1.0 0.0 0.0
outer loop
vertex $MAXX $MINY $MAXZ
vertex $MAXX $MAXY $MAXZ
vertex $MAXX $MAXY $MINZ
endloop
endfacet
endsolid WRAP
EOF
