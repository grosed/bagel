

rm example_2.o  particle.o  prior.o  dnorm.o theorem_2.o theorem_3.o theorem_1.o theorem_4.o bagel.o 
rm bagel

g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ example_2.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ particle.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ prior.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ dnorm.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ theorem_2.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ theorem_3.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ theorem_1.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ theorem_4.cpp
g++ -O3 -c -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ bagel.cpp

g++ -O3 -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I/home/grosed/DASS/ bagel_online.cpp example_2.o  particle.o  prior.o  theorem_1.o theorem_2.o theorem_3.o theorem_4.o dnorm.o bagel.o -o bagel
