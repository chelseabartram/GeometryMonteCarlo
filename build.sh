#!/bin/sh
clang++ -pthread -m64 -I/Users/cbartram/root-v5-34/include -m64 -L/Users/cbartram/root-v5-34/lib -c PSGenerator.cc
#-pthread -m64 -I/Users/cbartram/root-v5-34/include -m64 -L/Users/cbartram/root-v5-34/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl -rpath /Users/cbartram/root-v5-34/lib -lm -ldl -c PSGenerator.cc
clang++ -m64 -I/Users/cbartram/root-v5-34/include -m64 -L/Users/cbartram/root-v5-34/lib -c GeometryCalculator.cc
#-pthread -m64 -I/Users/cbartram/root-v5-34/include -m64 -L/Users/cbartram/root-v5-34/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl -rpath /Users/cbartram/root-v5-34/lib -lm -ldl -c GeometryCalculator.cc
clang++ -m64 -I/Users/cbartram/root-v5-34/include -m64 -L/Users/cbartram/root-v5-34/lib -c MonteCarlo.cc
#-pthread -m64 -I/Users/cbartram/root-v5-34/include -m64 -L/Users/cbartram/root-v5-34/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl -rpath /Users/cbartram/root-v5-34/lib -lm -ldl -c MonteCarlo.cc

clang++ -pthread -m64 -I/Users/cbartram/root-v5-34/include -m64 -L/Users/cbartram/root-v5-34/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl -rpath /Users/cbartram/root-v5-34/lib -lm -ldl -o PSTest MonteCarlo.o PSGenerator.o GeometryCalculator.o

