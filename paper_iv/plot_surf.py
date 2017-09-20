
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import scipy.io as sio

coord = sio.matlab.loadmat('./data/coord_power264.mat')['coord']

def draw_sphere(x, y, z, r):
    #draw sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = r * np.outer(np.cos(u), np.sin(v))+x
    y = r * np.outer(np.sin(u), np.sin(v))+y
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))+z
    return (x,y,z)

def plot_brain_surface(coord,size,color,spath):

    F=np.genfromtxt('./data/surface_face.csv',dtype='int')
    F=F-1
    V=np.genfromtxt('./data/surface_vertex.csv')
    FL = np.where(np.sum(V[F,0]>0,axis=1)==0)[0]
    FR = np.where(np.sum(V[F,0]>0,axis=1)==3)[0]

    #VL=np.genfromtxt('./data/surface_left_vertex_smooth.csv')
    #VR=np.genfromtxt('./data/surface_right_vertex_smooth.csv')
    #FL=np.genfromtxt('./data/surface_left_face_smooth.csv',dtype='int')
    #FL=FL-1
    #FR=np.genfromtxt('./data/surface_right_face_smooth.csv',dtype='int')
    #FR=FR-1

    #F = np.vstack([FL,FR])
    #V = np.vstack([VL,VR])

    Vmax = V.max()
    Vmin = V.min()




    fig=plt.figure(figsize=(20,10))
    ax1=plt.subplot2grid((2,4), (0, 0), rowspan=2,colspan=2,projection='3d')
    ax2=fig.add_subplot(243,projection='3d')
    ax3=fig.add_subplot(244,projection='3d')
    ax4=fig.add_subplot(247,projection='3d')
    ax5=fig.add_subplot(248,projection='3d')


    ax1.plot_trisurf(V[:,0],V[:,1],V[:,2],triangles=F,color=[0.8,0.8,0.8],linewidth=0, antialiased=True,edgecolor=[0.7,0.7,0.7],alpha=0.2)

    ax1.view_init(90,90)
    ax1.set_ylim(Vmin,Vmax)
    ax1.set_xlim(Vmin,Vmax)
    ax1.set_zlim(Vmin,Vmax)
    ax1.invert_yaxis()
    ax1.set_aspect('equal')
    ax1.axis('off')
    #ax1.set_position([0, 0, 1, 1])

    ax2.plot_trisurf(V[:,0],V[:,1],V[:,2],triangles=F[FR],color=[0.8,0.8,0.8],linewidth=0, antialiased=True,edgecolor=[0.7,0.7,0.7],alpha=0.2)
    ax2.view_init(0,0)
    ax2.set_ylim(Vmin,Vmax)
    ax2.set_xlim(Vmin,Vmax)
    ax2.set_zlim(Vmin,Vmax)
    #ax2.invert_xaxis()
    ax2.set_aspect('equal')
    ax2.axis('off')
    #ax2.set_position([0, 0, 1, 1])

    ax3.plot_trisurf(V[:,0],V[:,1],V[:,2],triangles=F[FL],color=[0.8,0.8,0.8],linewidth=0, antialiased=True,edgecolor=[0.7,0.7,0.7],alpha=0.2)
    ax3.view_init(0,180)
    ax3.set_ylim(Vmin,Vmax)
    ax3.set_xlim(Vmin,Vmax)
    ax3.set_zlim(Vmin,Vmax)
    ax3.set_aspect('equal')
    ax3.axis('off')
    #ax3.set_position([0,0,1,1])


    ax4.plot_trisurf(V[:,0],V[:,1],V[:,2],triangles=F[FR],color=[0.8,0.8,0.8],linewidth=0, antialiased=True,edgecolor=[0.7,0.7,0.7],alpha=0.2)
    ax4.view_init(0,180)
    ax4.set_ylim(Vmin,Vmax)
    ax4.set_xlim(0,Vmax-Vmin)
    ax4.set_zlim(Vmin,Vmax)
    ax4.set_aspect('equal')
    ax4.axis('off')
    #ax4.set_position([0, 0, 1, 1])


    ax5.plot_trisurf(V[:,0],V[:,1],V[:,2],triangles=F[FL],color=[0.8,0.8,0.8],linewidth=0, antialiased=True,edgecolor=[0.7,0.7,0.7],alpha=0.2)
    ax5.view_init(0,0)
    ax5.set_ylim(Vmin,Vmax)
    ax5.set_xlim(Vmin,Vmax)
    ax5.set_zlim(Vmin,Vmax)
    ax5.set_aspect('equal')
    ax5.axis('off')
    #ax5.set_position([0, 0, 1, 1])


    for n in range(0,coord.shape[0]):
        (x,y,z) = draw_sphere(coord[n,0],coord[n,1],coord[n,2],size[n])
        ax1.plot_surface(x, y, z, color=color[n],edgecolor='none')
        if coord[n,0]>0:
            ax2.plot_surface(x, y, z, color=color[n],edgecolor='none')
            ax4.plot_surface(x, y, z, color=color[n],edgecolor='none')
        elif coord[n,0]<0:
            ax3.plot_surface(x, y, z, color=color[n],edgecolor='none')
            ax5.plot_surface(x, y, z, color=color[n],edgecolor='none')

    plt.tight_layout()
    fig.savefig(spath,r='300')
