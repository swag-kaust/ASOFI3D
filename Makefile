# Default target to be used if not specified on the command line:
# that is, `make` is equivalent to `make all`.
.DEFAULT_GOAL := all

.PHONY : all
all :
	ls -l && cd src && $(MAKE) all

.PHONY : test
test :
	tests/test_01.sh
