.PHONY: clean anim gif octave build run correlation


CC  ?= gcc-15
CXX ?= g++-15

build:
	mkdir -p build
	cd build && cmake3 .. -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)
	$(MAKE) -C build

run:
	cd build/src && ./main 

anim: 
	cd build/src/script && python3 animation.py ./../output/test/

gif: 
	mkdir -p output
	cd build/src/script && python3 animation.py ./../output/test/ ./../../../output/anim.gif

octave: 
	cd build/src/script && octave --no-gui debug.m ./../output/test/sdf.csv ./../output/test/polyline.csv ./../output/test/sdf_params.csv

correlation: 
	cd build/src/script && python3 correlation.py ./../output/test/





clean:
	rm -rf build
