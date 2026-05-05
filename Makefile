.PHONY: clean anim gif octave build

BINARY := build/src/main
DONE   := build/src/output/test/.done


# prevent reruning when no changes
$(DONE): $(BINARY)
	cd build/src && ./main
	mkdir -p $(dir $(DONE))
	touch $(DONE)


anim gif octave:
	mkdir -p build
	cd build && cmake3 ..
	$(MAKE) -C build
	$(MAKE) _$@

_anim: $(DONE)
	cd build/src/script && python3 animation.py ./../output/test/
	cp build/src/output/test/sdf.json output/

_gif: $(DONE)
	mkdir -p output
	cd build/src/script && python3 animation.py ./../output/test/ ./../../../output/anim.gif

_octave: $(DONE)
	cd build/src/script && octave --no-gui debug.m ./../output/test/sdf.csv ./../output/test/polyline.csv



build:
	mkdir -p build
	cd build && cmake3 ..
	$(MAKE) -C build

clean:
	rm -rf build
