#! /bin/sh -
set -eu

module ()
{
    eval `/usr/local/Modules/bin/modulecmd bash $*`
}

# env for CCE build
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

touch compile_log.txt

# save module info
{
    module list
} |& tee -a compile_log.txt

# remove previous build, can cause build errors otherwise
{
    echo y | ./makenemo -n "${MY_CFG}" -r ORCA2_ICE_PISCES -m "${MY_ARCH}" -j 0 clean_config
} |& tee -a compile_log.txt

# Create directory for source code modifications
if [ ! -d "./MY_SRC" ]; then
    mkdir MY_SRC
else
    rm -rf ./MY_SRC
    mkdir MY_SRC
fi
cp -a "${PDAF_BIND_SRC}/pdaf_bindings/." ./MY_SRC/
cp -a "${PDAF_BIND_SRC}/nemo_src/." ./MY_SRC/
CWD=$(pwd)
export MY_SRC=${CWD}/MY_SRC

# create new NEMO config based on ORCA2_ICE_PISCES config,
# add source code modifications, build!
# WARNING: DO NOT use add_key option here as results in
# ~ factor of 10 increase in build time! Not sure why...
{
    ./makenemo -n "${MY_CFG}" -e "${MY_SRC}" -r ORCA2_ICE_PISCES -m "${MY_ARCH}" -j 32
} |& tee -a compile_log.txt
