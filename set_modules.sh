
echo '

From bash shell:

1. run the following commands to load the proper modules:
module purge
module load gcc/7.1.0-fasrc01 openmpi/3.1.1-fasrc01  cuda/10.0.130-fasrc01

2. set up the following variables:

export MOEBIUS_ROOT=_path_the_MAGIC_folder
export PYTHONPATH=$PYTHONPATH:$MOEBIUS_ROOT/BACKEND/SHOP:$MOEBIUS_ROOT/BACKEND/SCRIPTS
'



