#!/bin/sh

TEST_NAME_0="Number untouchables 10**4"
TEST_CORRECT_0=1212
./debug_main 10000 100; TEST_RESULT_0=$?
if [ $TEST_RESULT_0 -eq $TEST_CORRECT_0 ]; 
then
    echo "(PASS) $TEST_NAME_0"
else
    echo "(FAIL) $TEST_NAME_0, expected $TEST_CORRECT_0 but got $TEST_RESULT_0"
fi

TEST_NAME_1="memcheck 10**4"
TEST_CORRECT_1=1212
valgrind -q --error-exitcode=-1 ./debug_main 10000 100; TEST_RESULT_1=$?
if [ $TEST_RESULT_1 -eq $TEST_CORRECT_1 ]; 
then
    echo "(PASS) $TEST_NAME_1"
else
    echo "(FAIL) $TEST_NAME_1, expected $TEST_CORRECT_1 but got $TEST_RESULT_1"
fi