.PHONY: clean anim gif octave build run


run:
	cd build/src && ./main 

anim: 
	cd build/src/script && python3 animation.py ./../output/test/
	cp build/src/output/test/sdf.json output/

gif: 
	mkdir -p output
	cd build/src/script && python3 animation.py ./../output/test/ ./../../../output/anim.gif

octave: 
	cd build/src/script && octave --no-gui debug.m ./../output/test/sdf.csv ./../output/test/polyline.csv

correlation: 
	cd build/src/script && python3 correlation.py ./../output/test/



build:
	mkdir -p build
	cd build && cmake3 ..
	$(MAKE) -C build

clean:
	rm -rf build
