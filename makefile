# To include headers for python, for matplotlibcpp
run: main.cpp
	g++ -g -o run -O2 -std=c++11 -I/usr/include/python2.7 -lpython2.7 -L/usr/local/opt/openblas/lib -larmadillo -framework Accelerate  main.cpp
