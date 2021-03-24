#! /bin/sh -
set -eu

module ()
{
    eval `/usr/local/Modules/bin/modulecmd bash $*`
}

# env for CCE build
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

export PDAF_ARCH=cray_cce_mpi
touch compile_log.txt
cd src
make clean
make |& tee -a ../compile_log.txt
