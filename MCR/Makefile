all:
	@gcc \
		-Wall -Wextra \
		-Wpedantic \
		-O3 -funroll-loops \
		-g \
		rng.c \
		u512.S fp.S \
		mont.c \
		csidh.c \
		helper.c \
		main.c \
		-o main

bench:
	@gcc \
		-Wall -Wextra \
		-O3 -funroll-loops \
		-g \
		rng.c \
		u512.S fp.S \
		mont.c \
		csidh.c \
		helper.c \
		bench.c \
		-o main


debug:
	gcc \
		-Wall -Wextra \
		-g \
		rng.c \
		u512.S fp.S \
		mont.c \
		csidh.c \
		helper.c \
		main.c \
		-o main

clean:
	rm -f main

