#!/bin/bash
MYDIR=$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python $MYDIR/get_cline.pyc "$@"
