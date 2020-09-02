# Makefile for MacOS --note hardcoded path to XCode llvm-profdata
CFLAGS += -O3 -ffast-math -flto -DNDEBUG=1 $(PGOFLAGS)
LDLIBS = -lm 

PGO_GENFLAGS = -fprofile-generate --coverage
PGO_USEFLAGS = -fprofile-use
MUSCLE_PROF_FLAGS = -tree2 benchmark.nwk -cluster2 neighborjoining -distance1 kmer20_4 -diags -sv -maxiters 2

CPPSRC = $(sort $(wildcard *.cpp))
CPPOBJ	= $(subst .cpp,.o,$(CPPSRC))

default: muscle-pgo


$(CPPOBJ): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $@

muscle: $(CPPOBJ)
	$(CXX) $(CFLAGS) -o muscle $(CPPOBJ) $(LDLIBS)

muscle-pgo: 
	+$(MAKE) PGOFLAGS="$(PGO_GENFLAGS)" muscle
	./muscle -quiet -in benchmark.fa -out benchmark.faa $(MUSCLE_PROF_FLAGS)
	rm *.o muscle
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/llvm-profdata merge default*.profraw -o default.profdata
	+$(MAKE) PGOFLAGS="$(PGO_USEFLAGS)" muscle
