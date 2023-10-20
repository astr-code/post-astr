class input_var:

    def __init__(self,flowtype='not_defined',im=0,jm=0,km=0,ref_t=0.0,reynolds=0.0,mach=0.0,gridfile='not_defined'):
        self.flowtype = flowtype
        self.im = im
        self.jm = jm
        self.km = km
        self.ref_t = ref_t
        self.reynolds = reynolds
        self.mach = mach
        self.gridfile = gridfile

    def readvaule(self):

        input_file=read_inputfile_name()

        self.flowtype=read_flowtype(input_file)

        self.im,self.jm,self.km=read_dim(input_file)
        self.ref_t,self.reynolds,self.mach=read_flow_param(input_file)
        self.gridfile=read_gridfile(input_file)

        # return self.flowtype,self.im,self.jm,self.km

def read_inputfile_name():

    import sys
    print('command:',sys.argv)
    count=0
    for i in sys.argv:
        # print(sys.argv[count])
        if sys.argv[count] == '-input':
            filename=sys.argv[count+1]
            print('input file: ',filename)
            return filename
        count +=1

def read_flowfield_name():

    import sys
    print('command:',sys.argv)

    count=0

    for i in sys.argv:
        # print(sys.argv[count])
        if sys.argv[count] == '-ff':
            filename=sys.argv[count+1]
            print('flowfield file: ',filename)
            return filename
        count +=1

def read_sliceinfo():

    import sys
    print('command:',sys.argv)
    count=0
    for i in sys.argv:
        # print(sys.argv[count])
        if sys.argv[count] == '-slice':
            slicename=sys.argv[count+1]
            print('slice: ',slicename)
            return slicename
        count +=1

def read_view_var():

    import sys
    print('command:',sys.argv)
    count=0
    for i in sys.argv:
        # print(sys.argv[count])
        if sys.argv[count] == '-toview':
            varname=sys.argv[count+1:]
            print('variable to view: ',varname)
        count +=1

    return varname

def read_flowtype(filename):

    import re

    # filename = input("input file name: ")
    f = open(filename, "r")
    count=0
    for line_str in f:
      count +=1
      # line_str=f.readline()
      if count == 6:
        flowtype=re.sub('[\n]','',line_str)
        # print("flowtype: ",flowtype)
    f.close()

    return flowtype


def read_dim(filename):

    # filename = input("input file name: ")
    f = open(filename, "r")
    count=0
    for line_str in f:
      count +=1
      if count == 9:
        dim= [int(item) for item in line_str.split(',')]
        im=dim[0]
        jm=dim[1]
        km=dim[2]
    f.close()

    return im,jm,km

def read_flow_param(filename):

    import re

    f = open(filename, "r")

    count=0
    for line_str in f:
      count +=1
      # line_str=f.readline()
      if count == 24:
        item = re.sub('[d0\n]','',line_str)
        faram= [float(it) for it in item.split(',')]
        ref_t=faram[0]
        reynolds=faram[1]
        mach=faram[2]
        # print("    ref_t: ",ref_t)
        # print(" reynolds: ",reynolds)
        # print("     mach: ",mach)
      # print(count)
      # print(line_str)
      # print(line_str)
    f.close()

    return ref_t,reynolds,mach

def read_gridfile(filename):

    import re

    # filename = input("input file name: ")
    f = open(filename, "r")
    count=0
    for line_str in f:
      count +=1
      # line_str=f.readline()
      if count == 53:
        grifile=re.sub('[\n]','',line_str)
        # print("flowtype: ",flowtype)
    f.close()

    return grifile


def read_input_file():

    import re
    import sys
    print('command:',sys.argv)
    
    count=0
    for i in sys.argv:
        # print(sys.argv[count])
        if sys.argv[count] == '-input':
            filename=sys.argv[count+1]
            print('input file: ',filename)
        count +=1
    
    # filename = input("input file name: ")
    f = open(filename, "r")
    count=0
    for line_str in f:
      count +=1
      # line_str=f.readline()
      if count == 6:
        flowtype=line_str
        # print("flowtype: ",flowtype)
      if count == 9:
        dim= [int(item) for item in line_str.split(',')]
        im=dim[0]
        jm=dim[1]
        km=dim[2]
      if count == 24:
        item = re.sub('[d0\n]','',line_str)
        faram= [float(it) for it in item.split(',')]
        ref_t=faram[0]
        Reynolds=faram[1]
        Mach=faram[2]
        print("    ref_t: ",ref_t)
        print(" Reynolds: ",Reynolds)
        print("     Mach: ",Mach)
      # print(count)
      # print(line_str)
      # print(line_str)
    f.close()

    return flowtype,im,jm,km
