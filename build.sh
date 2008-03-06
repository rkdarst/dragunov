
for i in 01 02 03 10 11 12 13 ; do 

extra="-DFFid=$i"
echo compiling: "$extra"

gcc $* $extra -Wall -O2 -shared -fPIC -c dragunov_c.c
gcc -Wall -O2 -shared -fPIC -c -DMEXP=19937 -include SFMT-params.h SFMT.c

# With icc instead of gcc:
#icc -Wall -O2 -shared -fPIC -c dragunov_c.c
#icc -Wall -O2 -shared -fPIC -c -DMEXP=19937 -include SFMT-params.h SFMT.c

# link it together (gcc needed regardless of what compiler you use first)
gcc -Wall -O2 -shared -fPIC dragunov_c.o SFMT.o -o dragunov_${i}_c.so

done