all:
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99

slow: 
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99 -D slow

debug:
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99 -ggdb