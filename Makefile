# Default target to be used if not specified on the command line:
# that is, `make` is equivalent to `make all`.
.DEFAULT_GOAL := all

.PHONY : all
all :
	$(MAKE) --directory=src all

.PHONY : asofi3D
asofi3D :
	$(MAKE) --directory=src sofi3D

.PHONY : test
test :
	tests/check_test_env.sh
	tests/test_01.sh
	# tests/test_02.sh
	tests/test_03.sh
	tests/test_04.sh
	tests/test_05.sh
	tests/test_06.sh
	tests/test_07.sh
	tests/test_08.sh
