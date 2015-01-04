include base.mk

all: libsimplectic.so libsimplectic.a example

libsimplectic.so: simplectic.o
	$(CC) -shared -lm $(LDFLAGS) -o $@ $<

libsimplectic.a: simplectic.o
	$(AR) rcs libsimplectic.a simplectic.o

example: libsimplectic.a example.c
	$(CC) $(CFLAGS) -o $@ example.c -lm libsimplectic.a

clean:
	rm -f libsimplectic.a libsimplectic.so simplectic.o example
