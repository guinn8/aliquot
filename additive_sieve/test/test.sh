#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

$SCRIPT_DIR/../cli 10000 1000 > $SCRIPT_DIR/test.txt
cmp $SCRIPT_DIR/test.txt $SCRIPT_DIR/sn_10000.txt

if [ $? -eq 0 ]
then
    printf "\nTest Success! s(n) sieve is correct to 10000.\n\n"
    exit 0
else
    printf "\nTest Fail! s(n) sieve differs from known good list.\n\n"
    exit 1
fi