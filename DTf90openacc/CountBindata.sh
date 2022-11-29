#!/bin/bash


dir=bindata
if [ ! -d ${dir} ]; then
    mkdir ${dir}
fi

dir=output
if [ ! -d ${dir} ]; then
    mkdir ${dir}
fi

dir=figures
if [ ! -d ${dir} ]; then
    mkdir ${dir}
fi

dir=movies
if [ ! -d ${dir} ]; then
    mkdir ${dir}
fi

fstfile=`ls -1 bindata/unf*.dat 2>/dev/null | head -1`
echo $fstfile
declare -i fstnum=`echo  ${fstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`

lstfile=`ls -1 bindata/unf*.dat 2>/dev/null | tail -1`
echo $lstfile
declare -i lstnum=`echo  ${lstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`

echo ${fstnum} ${lstnum} > tmp
cat tmp | tr -d \\n  > control.dat
rm -fr tmp
