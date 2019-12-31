#!/bin/bash
lcov --capture --directory build/temp.linux-x86_64-3.7/src/gmpy2.gcno --output-file build/coverage.info
genhtml build/coverage.info --output-directory build/out
firefox build/out/index.html &
