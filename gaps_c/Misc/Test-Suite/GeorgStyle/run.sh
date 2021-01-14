#!/bin/ksh

set -u

NO_CONFIG_MF="foo"
export NO_CONFIG_MF

BASE="/stefan"

GAPC=$BASE/bin/gapc
GHC=ghc
MAKE=make
MAKEFLAGS=
PERL=perl

TEMP=./temp
GRAMMAR=../
LHS_DIR=..
RTLIB=$BASE/include/rtlib/
CPPFLAGS_EXTRA="-I../../../../ -O -DNDEBUG"
LDLIBS_EXTRA=""
RUN_CPP_FLAGS=""

KSH="ksh"


if [ -e $BASE/share/gapc/config.mf ]; then
  CONFIG_MF=$BASE/share/gapc/config.mf
else
  CONFIG_MF=$BASE/config/generic.mf
fi

err_count=0
succ_count=0
failed=0

FILTER=.

if [ $# == 2 ]; then
  FILTER=$2
fi

REF=../../../adpc-tng-logs/runs
if [ $# -ge 1 ]; then
  REF=$1
fi

mkdir -p $TEMP

if [ ! -d $TEMP/$REF ]; then
  echo Reference directory is no directory: $TEMP/$REF
  exit 1
fi

BASE_DIR="`pwd`"

. ./tool.sh

cd $BASE_DIR


mkdir -p $TEMP
cd $TEMP

echo include $CONFIG_MF > gapc_local.mf
printf RT_LDLIBS=\\n  >> gapc_local.mf
printf RT_LDLIBS04=\\n  >> gapc_local.mf

cmp_new_old_output()
{
  log1 s sort $1.$3.$4.out
  log1 a sed '/Answer/d' s

  log1 t sort $2/$1.$3.$4.out
  log1 b sed '/Answer/d' t

  log diff -u -w -B a b
}

check_new_old_eq()
{
  if [[ `echo $1$3$5 | grep $FILTER` != $1$3$5  ]]; then
    return
  fi

  # work around 1 sec timestamp filesystems ... WTF?!?
  sleep 1

  echo +------------------------------------------------------------------------------+
  failed=0
  temp=$failed

  cpp_base=${1%%.*}
  build_cpp $GRAMMAR/$1 $cpp_base $3 $2
  run_cpp $cpp_base $3 $4 $5
  cmp_new_old_output $cpp_base $REF $3 $5

  if [ $temp != $failed ]; then
    echo --++--FAIL--++--
    err_count=$((err_count+1))
  else
    echo OK
    succ_count=$((succ_count+1))
  fi
  echo +------------------------------------------------------------------------------+
}

run_check_feature()
{
  out=$1.$2.$3.$4
  echo ${*:5} $out
  "${@:5}" $out
  check_exit $?
}

check_feature()
{
  if [[ `echo $1$2 | grep $FILTER` != $1$2  ]]; then
    return
  fi
  echo +------------------------------------------------------------------------------+
  failed=0
  temp=$failed

  cpp_base=${1%%.*}
  build_cpp $GRAMMAR/$1 $cpp_base $2 ""
  run_cpp $cpp_base $2 $3 $4

  run_check_feature $cpp_base $2 $4 $5 "${@:6}"

  if [ $temp != $failed ]; then
    echo --++--FAIL--++--
    err_count=$((err_count+1))
  else
    echo OK
    succ_count=$((succ_count+1))
  fi
  echo +------------------------------------------------------------------------------+
}

check_compiler_output()
{
  if [[ `echo $1$2$3$4 | grep $FILTER` != $1$2$3$4  ]]; then
    return
  fi
  echo +------------------------------------------------------------------------------+

  cpp_base=${2%%.*}
  GRAMMAR=$1
  out=$cpp_base.$3.$4.gapc.log
  log_both $out  ${GAPC} ${GAPC_EXTRA} $GRAMMAR/$2 -o $cpp_base.cc -i $3

  failed=0
  temp=$failed

  echo ${*:5} $out
  "${@:5}" $out
  check_exit $?

  if [ $temp != $failed ]; then
    echo --++--FAIL--++--
    err_count=$((err_count+1))
  else
    echo OK
    succ_count=$((succ_count+1))
  fi
  echo +------------------------------------------------------------------------------+
}

check_feature_repeat_mean_var()
{
  if [[ `echo $1$4 | grep $FILTER` != $1$4  ]]; then
    return
  fi
  echo +------------------------------------------------------------------------------+
  failed=0
  temp=$failed

  cpp_base=${1%%.*}
  ###build_cpp $GRAMMAR/$1 $cpp_base $2


  out=$cpp_base.$2.$4.out
  for i in `../rand 100`; do
    GSL_RNG_SEED=$i ./$cpp_base -r 1000 $3 > log
    awk '/\[/ { sum++; array[$4]++; }
         END { for(i in array) print i, array[i]/(sum); } ' log | \
      sort -r -g -t' ' -k 2 | grep '^\[\] '
  done > z

#    /* printf("%s\t", $3); for (i=0; i<$1; i++) printf("#"); printf("\n"); */
  sort z | uniq -c | \
    awk '{ n+=$1; sum+=$1*$3; array[$1]=$3;
         }
         END { mean=sum/n; print "mean: ", mean;
               for (i in array) { a+=i*((array[i]-mean)^2); }
               print "sample variance: ", a/n; }' > l

  MEAN=`grep mean l | sed 's/^mean: //'`
  VAR=`grep sample l | sed 's/^[^:]\+: //'`

  echo Testing mean
  # 0.00059
  log ../../paraltest/fp_eq $MEAN 0.789261 0.004
  echo Testing variance 
  log ../../paraltest/fp_eq $VAR 2.87016e-05 0.00005


  if [ $temp != $failed ]; then
    echo --++--FAIL--++--
    err_count=$((err_count+1))
  else
    echo OK
    succ_count=$((succ_count+1))
  fi
  echo +------------------------------------------------------------------------------+
}

check_external()
{
  if [[ `echo $1 | grep $FILTER` != $1  ]]; then
    return
  fi
  echo +------------------------------------------------------------------------------+
  temp=$failed

  echo Testing $1
  log "${@:2}"

  if [ $temp != $failed ]; then
    echo --++--FAIL--++--
    err_count=$((err_count+1))
  else
    echo OK
    succ_count=$((succ_count+1))
  fi
  echo +------------------------------------------------------------------------------+
}

. ../config

. ../stats.sh

