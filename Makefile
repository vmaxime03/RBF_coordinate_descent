.PHONY: clean anim gif octave build build_release run run_release run_smoother correlation surface save colormap samples_error

CC   ?= gcc
CXX  ?= g++

DATA ?= run_output/test

CONFIG ?= smoother_config.json

build:
	mkdir -p build
	cd build && cmake .. -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)
	$(MAKE) -j8 -C build

build_release:
	mkdir -p build
	cd build && cmake .. -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_BUILD_TYPE=Release
	$(MAKE) -j8 -C build

run:
	cd build/src && ./main

run_smoother:
	cd build/src && ./smoother $(CURDIR)/$(CONFIG)


save:
	cp -r $(DATA) ./output/save_$(shell date +%Y%m%d_%H%M%S)

anim:
	cd build/src/script && python3 animation.py ./../../../$(DATA)/

gif:
	mkdir -p output
	cd build/src/script && python3 animation.py ./../../../$(DATA)/ ./../../../output/anim.gif

octave:
	cd build/src/script && octave debug_plot.m ./../../../$(DATA)/

surface:
	cd build/src/script && octave surface_plot.m ./../../../$(DATA)/

colormap:
	cd build/src/script && octave colormap_plot.m ./../../../$(DATA)/

samples_error:
	cd build/src/script && octave samples_plot.m ./../../../$(DATA)/

correlation:
	cd build/src/script && python3 correlation.py ./../../../$(DATA)/

clean:
	rm -rf build
