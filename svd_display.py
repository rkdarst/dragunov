import sys

import dragunov

S = dragunov.readUghfile(sys.argv[1])
S.makebox()
S.display()
