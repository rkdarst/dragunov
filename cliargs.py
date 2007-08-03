# Richard Darst, 2006.
# GPLv2+

"""This program is an easy command line parser.  The advantage is that
it is fast and easy, and that you don't have to specify what options
to use beforehand.

It returns a tuple consisting of:
- List of arguments
- Dictionary of options.
  - All options have the form '--xxx'
  - '--xxx' makes 'xxx' appear as a dict key with value of True.
  - '--xxx=yyy' makes 'xxx' appear in the dict with a vaule of 'yyy'.
  - '--' means 'no more options after this point'

See below for examples.

(rkd zgib.net)
"""

import sys

def get_cliargs(argv=None, multiopt=False, allowedargs=None):
    if argv is None:
        argv = sys.argv
    if type(argv) == str:
#        print "automatically converting to a list using split"
        argv = argv.split()
    argv = [ x for x in argv ]  # basically argv.copy(), if it existed
    options = { }
    args = [ ]
    while len(argv) > 0:
        nextarg = argv.pop(0)
        if nextarg == "--":  # standard for ending options
            args.extend(argv)
            break
        if nextarg[0:2] != "--":
            args.append(nextarg)
            continue
        nextarg = nextarg[2:]
        if nextarg.find("=") == -1:
            name = nextarg
            val = True
        else:
            name, val = nextarg.split("=", 1)
        if options.has_key(name):
            if not multiopt:
                raise Exception("already have this option: %s"%name)
            # handling of multiarg mode here
        if allowedargs and name not in allowedargs:
            print "Error-- arg %s found and not expected"%name
            break
        options[name] = val
    return args, options

if __name__ == "__main__":
    # Example use:
    sys.argv = ["program name", "--hi", "option", "--anoption", "--blah=xxx",
                # the "--" means "everything past this point is not an option"
                "--", "--notanoption"]
    args, options = get_cliargs()
    print args
    # --> ['program name', 'option', '--notanoption']
    print options
    # --> {'blah': 'xxx', 'hi': True, 'anoption': True}


    print
    print get_cliargs("--hi")
    # --> ([], {'hi': True})
    print get_cliargs("--hi --whathappen")
    # --> ([], {'hi': True, 'whathappen': True})
    print get_cliargs("--hi --withoption=one")
    # --> ([], {'withoption': 'one', 'hi': True})
    print get_cliargs("--hi option --anoption")
    # --> (['option'], {'hi': True, 'anoption': True})
    print get_cliargs("--hi option --anoption -- --notanoption")
    # --> (['option', '--notanoption'], {'hi': True, 'anoption': True})
    print get_cliargs(["--hi", "option", "--anoption", "--", "--notanoption"])
    # --> (['option', '--notanoption'], {'hi': True, 'anoption': True})
            
            
    
    
