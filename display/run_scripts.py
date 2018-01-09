import os


def main(args):
        
    script_path = 'display/graphs.r'
    file_path = args[1]

    command = "Rscript %s %s --vanilla"%(script_path, file_path)
    os.system(command)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))



