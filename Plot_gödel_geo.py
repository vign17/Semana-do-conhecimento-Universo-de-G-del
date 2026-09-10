import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

#Parametros do momento conjugado a t
k_0=1.5
k=[1,1.01,1.2,1.5]

#Cores de geodésicas
colors=['blue','yellow','orange','green']

#Limite de continuidade da solução para k0, caso luminal
lim=np.pi/(2*k_0)

#Limite de continuidade da solução para k, caso temporal
def lim2(En):
    limite=np.pi/(2*(En**2+1)**0.5)
    return limite

# Figura 3D
fig = plt.figure(figsize=(12, 10))

ax1 = fig.add_subplot(221, projection='3d')
ax2 = fig.add_subplot(222, projection='3d')
ax3 = fig.add_subplot(223, projection='3d')
ax4 = fig.add_subplot(224, projection='3d')


def plot_segmento_lum(limite,E,th_0,clr,ax):
    for it in range(1,3):
        tau = np.linspace(0, limite, 1000)
        
        # Curvas parametrizadas em tau, geodésicas luminais com mecanismo de extensão para segunda iteração, it=2
        r = np.arcsinh(np.sin(E*(tau+(it-1)*limite)))
        theta=th_0-np.arctan(2**0.5*np.tan(E*tau))-np.pi*(it-1)/2
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        t = -E*(tau+(it-1)*limite)+2**0.5*(np.arctan(2**0.5*np.tan(E*(tau-(it-1)*limite)))+np.pi*(it-1))

        ax.plot(x, y, t,color=clr)
        
def plot_segmento_temp(E,th_0,clr,ax):
    E2=E**2
    for it in range(1,3):
        tau = np.linspace(0, lim2(E), 10000)
        
        # Curvas parametrizadas em tau, geodésicas temporais com mecanismo de extensão para segunda iteração, it=2
        r = np.arcsinh(((E2-1)/(E2+1))**0.5*np.sin((E2+1)**0.5*(tau+(it-1)*lim2(E))))
        theta=th_0-np.arctan((2*E2/(E2+1))**0.5*np.tan((E2+1)**0.5*(tau)))-np.pi*(it-1)/2
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        t = -E*(tau+lim2(E)*(it-1))+2**0.5*(np.arctan((2*E2/(E2+1))**0.5*np.tan((E2+1)**0.5*(tau)-np.pi/2*(it-1)))+np.pi*(it-1))

        #Correção para prevenir artefatos no plot devido a erros numéricos
        jump = (
            (np.abs(np.diff(x)) > 0.1) |
            (np.abs(np.diff(y)) > 0.1) |
            (np.abs(np.diff(t)) > 0.1)
        )

        idx = np.where(jump)[0]

        x[idx + 1] = np.nan
        y[idx + 1] = np.nan
        t[idx + 1] = np.nan
        
        # Trajetória
        ax.plot(x, y, t,color=colors[clr])
        
        #Mecanismo de parada para caso trivial
        if(E==1):
            break

#Plot de geodésicas luminais para difentes ângulos iniciais, basta alterar parâmetro th_0 fornecido e range do for
def luminal(ax):    
    for i in range(1):
        plot_segmento_lum(lim,k_0,i*2*np.pi/5,'red',ax)
    
    # Eixos
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('t')

#Plot de geodésicas luminais para difentes ângulos iniciais, basta alterar parâmetro th_0 fornecido e range do for
def temporal(T,j,ax):    
    for i in range(1):
        plot_segmento_temp(T,i*2*np.pi/5,j,ax)

    # Eixos
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('t')

#Plot de geodésicas luminais para difentes ângulos iniciais, basta alterar parâmetro th_0 fornecido e range do for
for ax in [ax1, ax2, ax3, ax4]:
    ax.set_box_aspect((1,1,1))
    luminal(ax)

    for i in range(4):
        temporal(k[i], i,ax)

    legendas = [
        Line2D([0], [0], color='red', lw=2, label='Luminal'),
        *[
            Line2D([0], [0], color=colors[i], lw=2,
                label=fr'Parâmetro $E(0)/m = {k[i]}$')
            for i in range(4)
        ]
    ]
    
    #Determinação de limites de visualização do espaço
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    #Eixos
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('t')

# vista superior: XY
ax1.view_init(elev=90, azim=0)    

# vista frontal: XZ
ax2.view_init(elev=0, azim=90)    
# vistas oblíquas
ax3.view_init(elev=60, azim=270)     
ax4.view_init(elev=30, azim=45)

#Definição de legenda e título
fig.legend(handles=legendas,loc='lower right')
fig.suptitle(r'$\bf{Geodésicas \ radiais \ no \ universo \ de \ Gödel}$')

#Layout e reenderização
plt.tight_layout()
plt.show()