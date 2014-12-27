from tomopy.local.user_exception import IllegalArgumentError

def read_option(sys, help_string, nmin, nmax):
    if (len(sys.argv) < nmin or len(sys.argv) > nmax):
        raise IllegalArgumentError("number of arguments doesn't match.")
    else:
        option = sys.argv[1:]
        option_rdash = [s[2:] for s in option]
        option_dict = {}
        for opn in option_rdash:
            if opn == 'help':
                print(help_string)
                sys.exit()
            else:
                k, v = opn.split('=')
                option_dict[k] = v
    return option_dict
