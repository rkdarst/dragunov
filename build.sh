gcc -Wall -O2 -shared -fPIC -c dragunov_c.c
gcc -Wall -O2 -shared -fPIC dragunov_c.o -o dragunov_c.so

# With icc:
#icc -O2 -shared -fPIC -c dragunov_c.c
#gcc -O2 -shared -fPIC dragunov_c.o -o dragunov_c.so
