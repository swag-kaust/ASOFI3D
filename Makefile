.PHONY : all
all :
	$(MAKE) --directory=src all

.PHONY : asofi3D sofi3D
asofi3D sofi3D asofi3d sofi3d :
	$(MAKE) --directory=src sofi3D

.PHONY : clean
clean :
	$(MAKE) --directory=src clean

.PHONY : test
test :
	tests/check_test_env.sh
	tests/test_01.sh
	if [ -z "${CI}" ]; then tests/test_02.sh; fi
	tests/test_03.sh
	tests/test_04.sh
	tests/test_05.sh
	tests/test_06.sh
	tests/test_07.sh
	tests/test_08.sh
	tests/test_09.sh
	tests/test_10.sh
