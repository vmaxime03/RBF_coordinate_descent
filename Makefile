.PHONY: clean anim gif octave build run correlation surface save

CC  ?= gcc
CXX ?= g++

DATA ?= build/src/output/test

build:
	mkdir -p build
	cd build && cmake .. -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)
	$(MAKE) -C build

run:
	cd build/src && ./main

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

correlation:
	cd build/src/script && python3 correlation.py ./../../../$(DATA)/


clean:
	rm -rf build
