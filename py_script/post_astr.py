import h5io
import numpy as np

x=h5io.read_data('datin/grid.h5','x','kslice=50')
# print(type(x),x.ndim,x.shape)
y=h5io.read_data('datin/grid.h5','y','kslice=50')
# print(type(x),x.ndim,x.shape)
z=h5io.read_data('datin/grid.h5','z','kslice=0')
print(type(z),z.ndim,z.shape)

print(z)

# y=np.array([[1.0,3.0,5.0,7.0],[2.0,4.0,6.0,8.0]])
# print(type(y[0][0]),y.ndim,y.shape)
# print(y)



# import read_txt
# import h5io
# import primes
# inputdata=read_txt.input_var()
# inputdata.readvaule()
# # flowtype,im,jm,km = read_txt.read_input_file()

# print("flowtype: ",inputdata.flowtype,type(inputdata.flowtype))
# print("dimension: ",inputdata.im,'x',inputdata.jm,'x',inputdata.km,type(inputdata.km))
# print("ref_t: ",inputdata.ref_t,', reynolds: ',inputdata.reynolds,', mach: ',inputdata.mach,type(inputdata.mach))
# print("gridfile: ",inputdata.gridfile,type(inputdata.gridfile))

# print(h5io.h5_read1int.__doc__)
# print(primes.__doc__)
# print('all jobs done.')