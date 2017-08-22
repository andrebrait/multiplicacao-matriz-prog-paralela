all:
	gcc main.c -o bin/main

debug:
	gcc main.c -o bin/main -ggdb