#!/bin/bash

dirf=figures/
dirm=movies/

if [ $# -ne 1 ];then
    echo "usage:" $0 prefix
    exit 1
else
    prename=$1
fi

echo $prename

# first file
fstfile=`ls -1 ${dirf}/${prename}*.png  2>/dev/null | head -1`
echo $fstfile
declare -i fstnum=`echo  ${fstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`

ffmpeg -y -r 10 -start_number ${fstnum} -i ${dirf}${prename}%5d.png -b 6000k -vcodec wmv2 -pass 1 -r 10 -an ${dirm}ani${prename}.wmv

ffmpeg -y -r 10  -start_number ${fstnum} -i ${dirf}/${prename}%5d.png -vcodec libx264 -pix_fmt yuv420p -r 10 -an ${dirm}/ani${prename}.mp4

exit


