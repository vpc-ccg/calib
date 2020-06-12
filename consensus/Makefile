#!/usr/bin/make -rRf
my_mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
my_current_dir :=  $(patsubst %/,%,$(dir $(my_mkfile_path)))/

CXX?=g++
CMAKE?=cmake
cc_files?= calib_cons.cc
cc_flags?=
cc_args?= $(cc_flags)-std=c++11 -O3 -pthread -Ispoa_v1.1.3/include/ -Lspoa_v1.1.3/build/lib/ -lspoa
x:=$(shell $(CMAKE) --version | head -n1| cut -d " " -f3 | cut -c1-1)
$(info $x)
calib_cons: $(cc_files) spoa_v1.1.3/build/lib/libspoa.a 
	$(CXX) $(cc_files) $(cc_args) -o calib_cons

spoa_v1.1.3/build/lib/libspoa.a:
ifeq ($(x),3)
	@echo Using $(CMAKE)
	git clone --recursive https://github.com/rvaser/spoa spoa_v1.1.3
	cd spoa_v1.1.3; git checkout 1.1.3; mkdir build; cd build; $(CMAKE) -DCMAKE_BUILD_TYPE=Release ..;
	make -C spoa_v1.1.3/build
else
	$(error "cmake needs to be version 3")
endif
