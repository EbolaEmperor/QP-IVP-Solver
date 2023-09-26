JSONLIB=jsoncpp-1.9.5/build/lib/libjsoncpp.a
INC=-Ijsoncpp-1.9.5/include/json -Iinclude
LIB=-Ljsoncpp-1.9.5/build/lib -ljsoncpp -lquadmath
OPTIMIZE=-O2 -O3 -Ofast

solve: $(JSONLIB) solve.o IVP.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o functionFactory.o
	g++ solve.o IVP.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o functionFactory.o $(LIB) -o solve $(OPTIMIZE)

$(JSONLIB):
	cd jsoncpp-1.9.5 && mkdir build
	cd jsoncpp-1.9.5/build && cmake ..
	cd jsoncpp-1.9.5/build && make

solve.o: solve.cpp functions/*.h
	g++ -c solve.cpp $(INC) -Ifunctions $(OPTIMIZE)

functionFactory.o: src/functionFactory.cpp include/functionFactory.h
	g++ -c src/functionFactory.cpp $(INC) $(OPTIMIZE)

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp $(INC) $(OPTIMIZE)

IVP.o: src/IVP.cpp include/IVP.h src/matrix.cpp include/matrix.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/IVP.cpp $(INC) $(OPTIMIZE)

NonLinearSolver.o: src/NonLinearSolver.cpp include/NonLinearSolver.h src/matrix.cpp include/matrix.h
	g++ -c src/NonLinearSolver.cpp $(INC) $(OPTIMIZE)

RungeKutta.o: src/RungeKutta.cpp include/RungeKutta.h src/IVP.cpp include/IVP.h src/NonLinearSolver.cpp include/NonLinearSolver.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/RungeKutta.cpp $(INC) $(OPTIMIZE)

Polynomial.o: src/Polynomial.cpp include/Polynomial.h src/matrix.cpp include/matrix.h
	g++ -c src/Polynomial.cpp $(INC) $(OPTIMIZE)

clean:
	$(shell if [ -e solve ];then rm solve; fi)
	$(shell if [ -e result.txt ];then rm *.txt; fi)
	rm *.o
	rm -r jsoncpp-1.9.5/build
	cd report && make clean
	cd report && cd figures && mv 11-1.eps 11-1.epsback && rm *.eps 5-5.png && mv 11-1.epsback 11-1.eps
