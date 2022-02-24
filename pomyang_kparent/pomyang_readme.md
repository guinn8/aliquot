Pomerance-Yang Algorithm Readme
=============================

This project purpose it to implement the Pomerance-Yang algorithm efficiently, accurately, and scalable to high bounds. It includes a command line interface (cli) to interact with the algorithm dynamically and a test suite to help with accuracy. The documentation in the source files is extensive, pom_yang.c provides the best starting point.

Building and Running
--------------------

The cli can be build and ran with:

```bash
make cli && time ./cli --bound=$((10**8)) --seg_len=$((5 * $((10**5)))) --num_locks=$((10**7)) --preimage_count_bits=1 --num_thread=12
```

The tests will take maybe a couple a minutes to complete, successful if the program doesn't assert out. The tests can be ran with:

```bash
make test && ./test
```

Tested on Ubuntu 20.04.3 LTS.
