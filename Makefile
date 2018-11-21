.PHONY : all
all :
	$(MAKE) --directory=src all

.PHONY : asofi3D sofi3D asofi3d sofi3d
asofi3D sofi3D asofi3d sofi3d :
	$(MAKE) --directory=src sofi3D

.PHONY : clean
clean :
	$(MAKE) --directory=src $@

.PHONY : clobber
clobber :
	$(MAKE) --directory=src $@

.PHONY : test
test :
	tests/check_test_env.sh
	tests/test_01.sh
	if [ -n "${SLOW_TESTS}" ]; then tests/test_02.sh; fi
	tests/test_03.sh
	tests/test_04.sh
	tests/test_05.sh
	tests/test_06.sh
	tests/test_07.sh
	tests/test_08.sh
	tests/test_09.sh
	tests/test_10.sh

# Developer-level target, to check that one single translation unit
# compiles without any warnings from a compiler.
# Variable $(@F) contains the basename of the target %.o,
# so that the user could call, for example,
#     make src/snapmerge.o
# or
#     make snapmerge.o
# with the same result.
%.o:
	$(MAKE) --directory=src $(@F)
