all:
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99 -D IKJ

slow: 
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99 -D JKI

debug:
	gcc main.c dgemm-blocked.c -o bin/main -O1 -fstrict-aliasing -std=c99 -ggdb -D IKJ