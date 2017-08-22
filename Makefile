all:
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99

debug:
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99 -ggdb