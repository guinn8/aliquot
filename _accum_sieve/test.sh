#!/bin/sh
rm -f tmp.txt
valgrind ./sn_test 10000 100 > tmp.txt
if cmp tmp.txt data/sn_10000.txt; then
    echo "Successfully ennumerated s_n to 10000"
    # rm tmp.txt
else
    echo "Test failed"
fi

echo "time ./sn_test 10000000 100000 > /dev/null"
bash -c "time ./sn_test 10000000 100000 > /dev/null" 2>time.output  
cat time.output
rm time.output

