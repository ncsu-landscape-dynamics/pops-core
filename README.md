# SOD-modeling-cpp
recoding the model to create a c++ version of the SOD-model base on https://github.com/f-tonini/SOD-modeling.
This repository contains the c++ version scripts used to develop a stochastic landscape spread model of forest pathogen *P. ramorum*.

The reference paper: Ross K. Meentemeyer, Nik J. Cunniffe, Alex R. Cook, Joao A. N. Filipe, Richard D. Hunter, David M. Rizzo, and Christopher A. Gilligan 2011. Epidemiological modeling of invasion in heterogeneous landscapes: spread of sudden oak death in California (1990–2030). *Ecosphere* 2:art17. [http://dx.doi.org/10.1890/ES10-00192.1] (http://www.esajournals.org/doi/abs/10.1890/ES10-00192.1) 

## The files
The main.cpp contains the main program to run.

## To run the model(you can use Linux to run the model)
1.  Open an Linux shell and set the directory to the model folder.
2.  cd gdal-2.0.1
3.  sudo ./configure
4.  sudo make(The library will take about 1 hour to install)
5.  sudo make install
5.  cd ..
6.  sudo apt-get install libnetcdf-dev
7.  g++ -std=c++11 main.cpp Img.h Img.cpp Spore.h Spore.cpp -lgdal -lnetcdf_c++
8.  ./a.out > result.txt
