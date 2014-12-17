CXX=g++
CXXFLAGS=-pthread
NAME=spgrid_demo
all: Blocked_Copy_Helper.o Laplace_Helper.o PTHREAD_QUEUE.o SPGrid_Utilities.o main.o
	$(CXX) $(CXXFLAGS) -o $(NAME) \
		Blocked_Copy_Helper.o \
		Laplace_Helper.o \
		PTHREAD_QUEUE.o \
		SPGrid_Utilities.o \
		main.o
			
			
Blocked_Copy_Helper.o: Blocked_Copy_Helper.cpp
Laplace_Helper.o: Laplace_Helper.cpp
PTHREAD_QUEUE.o: PTHREAD_QUEUE.cpp
SPGrid_Utilities.o: SPGrid_Utilities.cpp
main.o: main.cpp
