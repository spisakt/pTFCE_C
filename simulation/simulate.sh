#!/bin/bash

# Determine RUNNING_SHELL; if SHELL is non-zero use that.
if [ -n "$SHELL" ]; then
    RUNNING_SHELL="$SHELL"
else
    if [ "$(uname)" = "Darwin" ]; then
        RUNNING_SHELL=/bin/bash
    else
        if [ -d /proc ] && [ -r /proc ] && [ -d /proc/$$ ] && [ -r /proc/$$ ] && [ -L /proc/$$/exe ] && [ -r /proc/$$/exe ]; then
            RUNNING_SHELL=$(readlink /proc/$$/exe)
        fi
        if [ -z "$RUNNING_SHELL" ] || [ ! -f "$RUNNING_SHELL" ]; then
            RUNNING_SHELL=$(ps -p $$ -o args= | sed 's|^-||')
            case "$RUNNING_SHELL" in
                */*)
                    ;;
                default)
                    RUNNING_SHELL=$(which "$RUNNING_SHELL")
                    ;;
            esac
        fi
    fi
fi

# Some final fallback locations
if [ -z "$RUNNING_SHELL" ] || [ ! -f "$RUNNING_SHELL" ]; then
    if [ -f /bin/bash ]; then
        RUNNING_SHELL=/bin/bash
    else
        if [ -f /bin/sh ]; then
            RUNNING_SHELL=/bin/sh
        fi
    fi
fi

if [ -z "$RUNNING_SHELL" ] || [ ! -f "$RUNNING_SHELL" ]; then
    printf 'Unable to determine your shell. Please set the SHELL env. var and re-run\\n' >&2
    exit 1
fi


THIS_DIR=$(DIRNAME=$(dirname "$0"); cd "$DIRNAME"; pwd)
THIS_FILE=$(basename "$0")
THIS_PATH="$THIS_DIR/$THIS_FILE"
SMOOTH=(0 0 0)
SNR=1
POSTSIGMA=1.5
TEST=0
USAGE="
usage: $0 [options]

Simulates data for pTFCE

-s S1 S2 S3  smoothing sigma for the three masks
-r SNR       signal to noise ratio
-p PS        smooting sigma for postprocessing
-t           dry run for testing
-h           print this help message and exit
"

#if which getopt > /dev/null 2>&1; then
#    OPTS=$(getopt bfhp:sut "$*" 2>/dev/null)
#    if [ ! $? ]; then
#        printf "%s\\n" "$USAGE"
#        exit 2
#    fi
#
#    eval set -- "$OPTS"
#
#    echo "${OPTS[@]}"

    while true; do
        case "$1" in
            -h)
                printf "%s\\n" "$USAGE"
                exit 2
                ;;
            -s)
                SMOOTH=($2 $3 $4)
                shift
                shift
                shift
                shift
                ;;
            -r)
                SNR=$2
                shift
                shift
                ;;
            -p)
                POSTSIGMA=$2
                shift
                shift
                ;;
            -t)
                TEST=1
                shift
                ;;
            --)
                shift
                break
                ;;
             *)
                if [[ $1 == '' ]]; then break; fi
                printf "ERROR: did not recognize option '%s', please try -h\\n" "$1"
                exit 1
                ;;
        esac
    done
#fi
#else
#    while getopts "bfhp:sut" x; do
#        case "$x" in
#            h)
#                printf "%s\\n" "$USAGE"
#                exit 2
#            ;;
#            b)
#                BATCH=1
#                ;;
#            f)
#                FORCE=1
#                ;;
#            p)
#                PREFIX="$OPTARG"
#                ;;
#            s)
#                SKIP_SCRIPTS=1
#                ;;
#            u)
#                FORCE=1
#                ;;
#            t)
#                TEST=1
#                ;;
#            ?)
#                printf "ERROR: did not recognize option '%s', please try -h\\n" "$x"
#                exit 1
#                ;;
#        esac
#    done
#fi

INDIR="${THIS_DIR}/noise"
OUTDIR="${THIS_DIR}/noise_${SMOOTH[0]}.${SMOOTH[1]}.${SMOOTH[2]}_${SNR}_${POSTSIGMA}"
echo "$OUTDIR"

if [ $TEST == 1 ]; then
  echo "${SMOOTH[@]}"
  echo "$SNR"
  echo "$POSTSIGMA"
  echo "testrun"
  #exit 1
fi

if [ $TEST == 0 ]; then mkdir -p $OUTDIR; fi
FILES=($(ls ${THIS_DIR}/noise))

for f in ${FILES[@]}
do
    echo "${f}"

    if [ $TEST == 1 ]; then
    echo "test"
    exit 1
    fi

    fslmaths ${INDIR}/${f} -mul ${THIS_DIR}/sim_mask_150.nii.gz ${OUTDIR}/n0.nii.gz
    fslmaths ${INDIR}/${f} -s ${SMOOTH[0]} -mul ${THIS_DIR}/sim_mask_90.nii.gz ${OUTDIR}/n1.nii.gz
    fslmaths ${OUTDIR}/n1.nii.gz -div `fslstats ${OUTDIR}/n1.nii.gz -S` ${OUTDIR}/n1.nii.gz
    fslmaths ${INDIR}/${f} -s ${SMOOTH[1]} -mul ${THIS_DIR}/sim_mask_60.nii.gz ${OUTDIR}/n2.nii.gz
    fslmaths ${OUTDIR}/n2.nii.gz -div `fslstats ${OUTDIR}/n2.nii.gz -S` ${OUTDIR}/n2.nii.gz
    fslmaths ${INDIR}/${f} -s ${SMOOTH[2]} -mul ${THIS_DIR}/sim_mask_30.nii.gz ${OUTDIR}/n3.nii.gz
    fslmaths ${OUTDIR}/n3.nii.gz -div `fslstats ${OUTDIR}/n3.nii.gz -S` ${OUTDIR}/n3.nii.gz

    fslmaths ${OUTDIR}/n0.nii.gz -add ${OUTDIR}/n1.nii.gz -add ${OUTDIR}/n2.nii.gz -add ${OUTDIR}/n3.nii.gz ${OUTDIR}/nx.nii.gz
    fslmaths ${THIS_DIR}/sim_signal.nii.gz -mul ${SNR} -add ${OUTDIR}/nx.nii.gz ${OUTDIR}/nxs.nii.gz

    fslmaths ${OUTDIR}/nx.nii.gz -s ${POSTSIGMA} ${OUTDIR}/n.nii.gz
    fslmaths ${OUTDIR}/n.nii.gz -div `fslstats ${OUTDIR}/n.nii.gz -S` ${OUTDIR}/no.nii.gz
    fslmaths ${THIS_DIR}/sim_signal.nii.gz -mul ${SNR} -add ${OUTDIR}/no.nii.gz ${OUTDIR}/ns.nii.gz
    fslmaths ${OUTDIR}/nxs.nii.gz -s ${POSTSIGMA} ${OUTDIR}/ns2.nii.gz
    fslmaths ${OUTDIR}/ns2.nii.gz -div `fslstats ${OUTDIR}/n.nii.gz -S` ${OUTDIR}/ns2.nii.gz

    fslroi ${OUTDIR}/ns.nii.gz ${OUTDIR}/s_${f} 30 90 30 90 30 90
    fslroi ${OUTDIR}/ns2.nii.gz ${OUTDIR}/sm_s_${f} 30 90 30 90 30 90
    fslroi ${OUTDIR}/no.nii.gz ${OUTDIR}/${f} 30 90 30 90 30 90
done

rm ${OUTDIR}/n0.nii.gz ${OUTDIR}/n1.nii.gz ${OUTDIR}/n2.nii.gz ${OUTDIR}/n3.nii.gz ${OUTDIR}/n.nii.gz ${OUTDIR}/ns.nii.gz ${OUTDIR}/nx.nii.gz ${OUTDIR}/ns2.nii.gz ${OUTDIR}/nxs.nii.gz

``
