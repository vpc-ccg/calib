#!/usr/bin/make -rRf
my_mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
my_current_dir :=  $(patsubst %/,%,$(dir $(my_mkfile_path)))/

cc?=g++
cc_files?= calib_cons.cc
cc_flags?=
cc_args?= $(cc_flags)-std=c++11 -O3 -pthread -Ispoa_v1.1.3/include/ -Lspoa_v1.1.3/build/lib/ -lspoa



CMAKE=cmake
ifeq ($(shell $(CMAKE) --version | awk '{split($NF, a, "."); if (int(a[1]) == 3 && int(a[2]) >= 2 || int(nums[1]) > 3) {print 100}}' 2>&1),100)
	@echo $(shell which $(CMAKE)) is fine
	CMAKE=cmake
else
	# @echo $(shell which $(CMAKE)) is not fine
	CMAKE=$(my_current_dir)cmake-3.12.0/bin/cmake
endif

calib_cons: $(cc_files)
	$(cc) $(cc_files) $(cc_args) -o calib_cons


spoa_v1.1.3/build/: $(CMAKE)
	@echo Using $(CMAKE)
	git clone --recursive https://github.com/rvaser/spoa spoa_v1.1.3
	cd spoa_v1.1.3; git checkout 1.1.3; mkdir build; cd build; $(CMAKE) -DCMAKE_BUILD_TYPE=Release ..;
	make -C spoa_v1.1.3/build

$(my_current_dir)cmake-3.12.0/bin/cmake:
	@echo Making cmake_3_12_0
	wget https://cmake.org/files/v3.12/cmake-3.12.0.tar.gz
	tar -xzvf cmake-3.12.0.tar.gz
	rm cmake-3.12.0.tar.gz
	cd cmake-3.12.0 && ./bootstrap;
	make -C cmake-3.12.0
