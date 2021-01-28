.PHONY: guard all fit fit-release fit-debug zeroer draw_graph

lib-src := $(wildcard src/*.cpp)
minuit-flags := ./minuit/lib/libminuit.a -lgfortran
#opt-flags := -O3
opt-flags := -Og -ggdb
flags := $(opt-flags) -Wall -Wextra -Werror -Wno-sign-compare -Wno-unused

guard:
	@echo "No-no-no. Write command explicitly please."

all: fit draw_graph

fit: fit-release

fit-release: Release
	@$(MAKE) -C Release install

fit-debug: Debug
	@$(MAKE) -C Debug install

Release:
	@mkdir $@
	@cmake -S . -B $@ -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=$@

Debug:
	@mkdir $@
	@cmake -S . -B $@ -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=$@

zeroer:
	g++ -I include zeroer.cpp -std=gnu++17 $(flags) $(lib-src) $(minuit-flags) -o $@

zerofit: zerofit.cpp
	g++ $^ `root-config --cflags --libs` -O3 -o $@

show-timesym-onemeasure:
	@./fit -i parms/timesym-onemeasure.parm -mf1 -p0,3 -p1,4 -p2,5 -R1e8 -S1e15 -I <<<"end"

show-spatsym-timesym-onemeasure:
	@./fit -i parms/spatsym-timesym-onemeasure.parm -mf1 --matrix-def spatsym -p0,3 -p1,4 -p2,5 -R1e12 -S1e20 -I <<<"end"

show-spatsym-timesym-manymeasures:
	@./fit -i parms/spatsym-timesym-manymeasures.parm -mf1 -mf2 -mfull -mfull -mfull -mmult -mfull -mfull -mmult -p0,3 -p1,4 -p2,5 -R1e6 -S1e12 --matrix-def spatsym -I <<<"end"

show-test-model:
	@./fit -i parms/test-model.parm -mf1 -mf2 -mfull -mfull -mfull -mmult -p0,3 -p1,4 -p2,5 -R1e16 -S1e-3 --all-cells --matrix-def spatsym -I <<<"end"

show-timesym-manymeasures:
	@./fit -i parms/timesym-manymeasures.parm -mf1 -mf2 -mfull -mfull -mfull -mfull -mmult -p0,3 -p1,4 -p2,5 -R1 -S1e10 -I <<<"end"

show-timesym-manymeasures-allpairs:
	@./fit -i parms/timesym-manymeasures-allpairs.parm -mf1 -mf2 -mfull -mfull -mfull -mmult -mfull -mfull -mfull -mfull -mfull -p0,3 -p1,4 -p2,5 -R1e7 -S1e15 --all-cells -I <<<"end"

show-manymeasures:
	@./fit -i parms/manymeasures.parm -mf1 -mf2 -mfull -mfull -mfull -mfull -mmult -p0,3 -p1,4 -p2,5 -R 1e15 -S1e-3 -I <<<"end"

show-the-fall-of-ghosts:
	@./fit -N3 --matrix-def spatsym -mf1 -mf2 -mfull -mfull -mfull -mmult --all-cells -p0,3 -p1,4 -p2,5 --load-matrix ~/bosons-with-symmetry-N3-the-fall-of-ghosts.txt -I <<<"end"

show-baseline:
	@./fit -N3 --matrix-def spatsym -mf1 -mf2 -mfull -mfull -mfull -mmult --all-cells -p0,3 -p1,4 -p2,5 --load-matrix ~/boltzman-with-symmetry-N3-baseline.txt -I <<<"end"
