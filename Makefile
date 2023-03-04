export CXXFLAGS = -Wall -Wextra
export LDFLAGS =

debug: CXXFLAGS += -O2 -fsanitize=address -fno-omit-frame-pointer
debug: LDFLAGS += -fsanitize=address

release:
	rm -rf $@
	cmake -S . -B $@ -G Ninja -DCMAKE_BUILD_TYPE=Release
	cmake --build $@

debug:
	rm -rf $@
	cmake -S . -B $@ -G Ninja -DCMAKE_BUILD_TYPE=Debug
	cmake --build $@

clean:
	rm -rf release debug

.PHONY: release debug clean
