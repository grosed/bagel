

rm example_2.o  particle.o  prior.o  theorem_2.o  theorem_3.o
rm praxi

g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ example_2.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ particle.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ prior.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ dnorm.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ theorem_2.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ theorem_3.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ theorem_1.cpp


g++ -O3 -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ praxi.cpp example_2.o  particle.o  prior.o  theorem_2.o theorem_3.o dnorm.o -o praxi

g++ -O3 -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ bagel_online.cpp example_2.o  particle.o  prior.o  theorem_2.o theorem_3.o dnorm.o -o bagel
