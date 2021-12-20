#!/bin/sh

BOUND=100000

ACCUM_FILE=sn_accum_$BOUND.txt
SIEVE_FILE=sn_sieve_$BOUND.txt

./sn_accum $BOUND 10 > $ACCUM_FILE
./sn_sieve $BOUND > $SIEVE_FILE

if cmp $ACCUM_FILE $SIEVE_FILE; then
    echo "Successfully ennumerated s_n to $BOUND"
else
    echo "Test failed"
fi

# rm $ACCUM_FILE
# rm $SIEVE_FILE