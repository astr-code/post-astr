import read_txt
import numpy as np

def tecxy(flowfield,slice_cut,tecname):

    import h5io
    from tecio import tecio
    import gradsolver as gs
    
    # slice_cut=read_txt.read_sliceinfo()
    # var2view =read_txt.read_view_var()
    # flowfield=read_txt.read_flowfield_name()
    
    # print(var2view,type(var2view))
    
    x=h5io.read_data('datin/grid.h5','x',slice_cut)
    print(type(x),type(x[0][0]),x.ndim,x.shape)
    y=h5io.read_data('datin/grid.h5','y',slice_cut)
    
    # ro=h5io.read_data(flowfield,'ro',slice_cut)
    # u1=h5io.read_data(flowfield,'u1',slice_cut)
    # u2=h5io.read_data(flowfield,'u2',slice_cut)
    # u3=h5io.read_data(flowfield,'u3',slice_cut)
    # p =h5io.read_data(flowfield,'p',slice_cut)
    # t =h5io.read_data(flowfield,'t',slice_cut)
    # sp01 =h5io.read_data(flowfield,'sp01',slice_cut)
    # sp02 =h5io.read_data(flowfield,'sp02',slice_cut)
    
    # omegaz=gs.omegaz(u1,u2,x,y)
    
    # du=gs.df2d(u1,x,y)
    # dv=gs.df2d(u2,x,y)
    
    # omegaz=dv[0,:,:]-du[1,:,:]
    
    # print(type(omegaz),omegaz.ndim,omegaz.shape)
    sigmaxx=h5io.read_data(flowfield,'sigmaxx',slice_cut)
    sigmaxy=h5io.read_data(flowfield,'sigmaxy',slice_cut)
    sigmaxz=h5io.read_data(flowfield,'sigmaxz',slice_cut)
    sigmayy=h5io.read_data(flowfield,'sigmayy',slice_cut)
    sigmayz=h5io.read_data(flowfield,'sigmayz',slice_cut)
    sigmazz=h5io.read_data(flowfield,'sigmazz',slice_cut)

    tecio.w2d8var(tecname,x,'x',y,'y',sigmaxx,'sigmaxx',sigmaxy,'sigmaxy',sigmaxz,'sigmaxz',sigmayy,'sigmayy',sigmayz,'sigmayz',sigmazz,'sigmazz')
    
    # tecio.w2d5var(tecname,x,'x',y,'y',omegaz,'omegaz',sp01,'Y1',sp02,'Y2')


slice_cut=read_txt.read_sliceinfo()
filename=read_txt.read_flowfield_name()
tecname='tecxy.plt'

tecxy(filename,slice_cut,tecname)

# for i in range(1,21):

#     filename='outdat/flowfield%0.4d.h5'%i
#     tecname='tecxy%0.4d.plt'%i
#     # print(filename,tecname)

#     tecxy(filename,slice_cut,tecname)


# count=0
# for varname in var2view:

#     print('data to view:', var2view[count],type(var2view[count]))

#     data[count]=h5io.read_data(flowfield,varname,slice_cut)

#     print(type(data[count]),data[count].ndim,data[count].shape)

#     count +=1


# print(type(x),x.ndim,x.shape)
# z=h5io.read_data('datin/grid.h5','z',slice_cut)
# print(type(z),z.ndim,z.shape)


# tecio.w2d2var('test1.plt',x,'x',y,'y')
# tecio.w2d3var('test2.plt',x,'x',y,'y',z,'z')


# print(z)
