#!/bin/sh
for i in `seq 9 12`
do
    BOUND=$(echo 10^$i | bc)
    # for j in `seq 3 $i`
    # do
        BLOCK=$(echo 10^5 | bc)
        CMD="./sn $BOUND $BLOCK"
        echo $CMD
        nice -n-20 /usr/bin/time -f%C,%e -o timings.csv -a $CMD > /dev/null 
    # done
   
done

