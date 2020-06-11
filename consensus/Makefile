#!/usr/bin/make -rRf
my_mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
my_current_dir :=  $(patsubst %/,%,$(dir $(my_mkfile_path)))/

CXX?=g++
cc_files?= calib_cons.cc
cc_flags?=
cc_args?= $(cc_flags)-std=c++11 -O3 -pthread -Ispoa/include/ -Lspoa/build/lib/ -lspoa

CMAKE?=cmake

calib_cons: $(cc_files) spoa/build/lib/libspoa.a
	$(CXX) $(cc_files) $(cc_args) -o calib_cons

spoa/build/lib/libspoa.a:
	mkdir -p spoa/build
	cd spoa/build && $(CMAKE) -DCMAKE_BUILD_TYPE=Release ..
	$(MAKE) -C spoa/build
