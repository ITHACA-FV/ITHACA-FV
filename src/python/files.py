import os
import re
import shutil
from tempfile import mkstemp

def sed(what, pattern, replace, source, dest=None):
    """Reads a source file and writes the destination file.

    In each line, replaces pattern with replace.

    Args:
    	what (str): pattern of the line to be replaced
        pattern (str): pattern to match (can be re.pattern)
        replace (str): replacement str
        source  (str): input filename
        dest (str):   destination filename, if not given, source will be over written.        
    """

    fin = open(source, 'r')

    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')

    for line in fin:
        if what in line:
            out = re.sub(pattern, replace, line)
            fout.write(out)
        else:
            fout.write(line)

    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E

    fin.close()
    fout.close()
    fout.close()

    if not dest:
        shutil.move(name, source) 

def sed_variable(what, source, new_value, dest=None):
    """Reads a source file and writes the destination file.

    In each line, replaces pattern with replace.

    Args:
        what (str): pattern of the line to be replaced
        source  (str): input filename
        new_value (str): replacement str
        dest (str):   destination filename, if not given, source will be over written.        
    """

    fin = open(source, 'r')

    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')

    for line in fin:
        if what in line:
            a = line.split();
            a[1] = str(new_value)+";\n"
            l = " ".join(str(x) for x in a)
            fout.write(l)
        else:
            fout.write(line)

    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E

    fin.close()
    fout.close()
    fout.close()

    if not dest:
        shutil.move(name, source)

def read_variable(what, source):
    """Reads a source file and returns the value of a variable.

    Args:
        what (str): pattern of the line to be returned
        source  (str): input filename.       
    """

    fin = open(source, 'r')
    fd, name = mkstemp()

    for line in fin:
        if what in line:
            a = line.split();
            return(int(a[1][:-1]))

    fin.close()

    if not dest:
        shutil.move(name, source) 