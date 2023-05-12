problem2a: problem2a.o problem.o map.o stack.o pq.o
	gcc -Wall -o problem2a problem2a.o problem.o map.o stack.o pq.o -g -lm

problem2b: problem2b.o problem.o map.o stack.o pq.o
	gcc -Wall -o problem2b problem2b.o problem.o map.o stack.o pq.o -g -lm

problem2d: problem2d.o problem.o map.o stack.o pq.o
	gcc -Wall -o problem2d problem2d.o problem.o map.o stack.o pq.o -g -lm

problem2e: problem2e.o problem.o map.o stack.o pq.o
	gcc -Wall -o problem2e problem2e.o problem.o map.o stack.o pq.o -g -lm
	

problem2a.o: problem2a.c problem.h map.h stack.h pq.h
	gcc -c problem2a.c -Wall -g -lm

problem2b.o: problem2b.c problem.h map.h stack.h pq.h
	gcc -c problem2b.c -Wall -g -lm

problem2d.o: problem2d.c problem.h map.h stack.h pq.h
	gcc -c problem2d.c -Wall -g -lm

problem2e.o: problem2e.c problem.h map.h stack.h pq.h
	gcc -c problem2e.c -Wall -g

problem.o: problem.c problem.h map.h stack.h pq.h
	gcc -c problem.c -Wall -g -lm

map.o: map.c map.h stack.h pq.h
	gcc -c map.c -Wall -g -lm

stack.o: stack.c stack.h
	gcc -c stack.c -Wall -g -lm

pq.o: pq.c pq.h
	gcc -c pq.c -Wall -g -lm
