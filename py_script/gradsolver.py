import numpy as np
from fgrad import fgradsolver

def df2d(f,x,y):

    im=len(x)
    jm=len(x[0])

    print('dimension of the array f:',im,jm)

    # df=fgradsolver.grad_xy(im,jm,f,x,y)
    df=fgradsolver.grad_xy(im,jm,f,x,y)

    # print(type(df),df.ndim,df.shape)

    return df

def omegaz(u,v,x,y):

    im=len(u)
    jm=len(u[0])


    du=fgradsolver.grad_xy(im,jm,u,x,y)  
    dv=fgradsolver.grad_xy(im,jm,v,x,y)

    omega=dv[0,:,:]-du[1,:,:]

    return omega

    # df=np.gradient(f, dx, dy)

    # print(range(len(x)),len(x))
    # print(x[0][0],x[260][0],x[0][740])
    # for i in range(len(x)):
    # dx=np.gradient(x)
    # dy=np.gradient(y)

    # dx1=np.arange(x.size)
    # dx2=np.arange(x.size)
    # dy1=np.arange(y.size)
    # dy2=np.arange(y.size)
    
    # dx1=dx[0]
    # dx2=dx[1]
    # dy1=dy[0]
    # dy2=dy[1]

    # dx1=dx[0]
    # dx2=dx[1]
    # dy1=dy[0]
    # dy2=dy[1]


    # print(dx2[0][0],type(dx1),dx1.ndim,dx1.shape)

    # for i in range(len(x)):
    #     for j in range(len(x[0])):
    #         var1=dx2[i][j]*dy1[i][j]-dx1[i][j]*dy2[i][j]
    # print(len(x),len(dx1))
    # print(dx1)
    # print(x)
    # dx1=
    # print(ff)
    # print(dff)
    # # =np.gradient(x)
    # print(type(ff),type(dff[0]))
    # return df