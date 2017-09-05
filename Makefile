MODULE_TOPDIR = ../..

PGM = sod-cpp

LIBES = $(RASTERLIB) $(GISLIB) $(MATHLIB) $(VECTORLIB)
DEPENDENCIES = $(RASTERDEP) $(GISDEP) $(VECTORDEP)
# $(NETCDFLIBS) is only C
EXTRA_LIBS = $(GDALLIBS) -lnetcdf_c++ $(OMPLIB)
EXTRA_CFLAGS = $(GDALCFLAGS) -std=c++11 -Wall -Wextra -fpermissive $(OMPCFLAGS) $(VECT_CFLAGS)
EXTRA_INC = $(VECT_INC)

include $(MODULE_TOPDIR)/include/Make/Module.make

LINK = $(CXX)

ifneq ($(strip $(CXX)),)
default: cmd
endif

# to compile as standalone, comment out the lines above
# TODO: add also -Wshadow -Wsign-conversion
standalone:
	g++ -std=c++11 -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused main.cpp Img.h Img.cpp Spore.h Spore.cpp -lgdal -lnetcdf_c++
