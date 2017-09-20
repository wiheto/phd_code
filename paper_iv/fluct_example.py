
    import numpy as np
    import teneto
    import matplotlib.pyplot as plt
    plt.rcParams['image.cmap'] = 'gist_gray'

    #Brute force make a 3D matrix, 12 time points long with repeating patterns
    A=np.zeros((4,4,1))
    B=np.zeros((4,4,1))
    A[0,1,0]=1
    A[2,3,0]=1
    B[0,3,0]=1
    B[0,1,0]=1

    A_nonflex=np.concatenate([A,B,np.tile(A,5),np.tile(B,2),np.tile(A,3)],2)

    #Brute force make a 3D matrix, 12 time points long with lots of different edges
    A=np.zeros((4,4,1))
    B=np.zeros((4,4,1))
    C=np.zeros((4,4,1))
    D=np.zeros((4,4,1))
    A[0,1,0]=1
    A[2,3,0]=1
    B[0,3,0]=1
    B[0,1,0]=1
    C[1,3,0]=1
    C[1,2,0]=1
    D[0,1,0]=1
    D[0,2,0]=1
    A_vol=np.concatenate([A,B,B,C,D,A,D,C,B,B,D,C],2)


    A_var=np.concatenate([A,A,D,D,D,B,B,B,B,C,C,C],2)


    fig,ax = plt.subplots(3,1)


    ax[0] = teneto.plot.slice_plot(A_nonflex,ax[0],vlabs=range(1,5),dlabs=list(map(str,range(1,13))))
    ax[0].set_title('A',loc='left')
    ax[0].set_ylabel('nodes')

    ax[1] = teneto.plot.slice_plot(A_vol,ax[1],vlabs=range(1,5),dlabs=list(map(str,range(1,13))))
    ax[1].set_title('B',loc='left')
    ax[1].set_ylabel('nodes')

    ax[2] = teneto.plot.slice_plot(A_var,ax[2],vlabs=range(1,5),dlabs=list(map(str,range(1,13))))
    ax[2].set_xlabel('time')
    ax[2].set_ylabel('nodes')
    ax[2].set_title('C',loc='left')
    fig.show()
    fig.tight_layout()

    fig.savefig('./examples/figures/fluctvol_example.pdf')
