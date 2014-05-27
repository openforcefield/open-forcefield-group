import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as MT
import matplotlib.lines as L
import matplotlib.cm as CM
import matplotlib.colors as C
import matplotlib.patches as PA
import pandas as pd
import mdtraj as md
from dipeptide_parameters import *

reference = pd.read_csv("./experimental_data/baldwin_table1_populations.csv", index_col=0)

data = []
for (ff, water, seq) in products:
    try:
        aa = seq.split("_")[1]
        t = md.load("./dcd/%s_%s_%s.dcd" % (ff, water, seq), top="./pdbs/%s.pdb" % (seq))
    except:
        continue
    phi = md.compute_phi(t)[1][:, 0] * 180 / np.pi
    psi = md.compute_psi(t)[1][:, 0] * 180 / np.pi
    ass = assign(phi, psi)
    populations = pd.Series({"PPII":0.0, "beta":0.0, "alpha":0.0, "other":0.0})
    populations += ass.value_counts(normalize=True)
    data.append([ff, water, aa, populations["PPII"], populations["beta"], populations["alpha"]])

data = pd.DataFrame(data, columns=["ff", "water", "aa", "PPII", "beta", "alpha"])


def plotSimplex(points, fig=None, vertexlabels=['PPII','beta','alpha'], **kwargs):
    """
    Plot Nx3 points array on the 3-simplex 
    (with optionally labeled vertices) 
    
    kwargs will be passed along directly to matplotlib.pyplot.scatter    

    Returns Figure, caller must .show()
    """
    if(fig == None):        
        fig = plt.figure()
    # Draw the triangle
    l1 = L.Line2D([0, 0.5, 1.0, 0], # xcoords
                  [0, np.sqrt(3) / 2, 0, 0], # ycoords
                  color='k')
    fig.gca().add_line(l1)
    fig.gca().xaxis.set_major_locator(MT.NullLocator())
    fig.gca().yaxis.set_major_locator(MT.NullLocator())
    # Draw vertex labels
    fig.gca().text(-0.05, -0.05, vertexlabels[0])
    fig.gca().text(1.05, -0.05, vertexlabels[1])
    fig.gca().text(0.5, np.sqrt(3) / 2 + 0.05, vertexlabels[2])
    # Project and draw the actual points
    projected = projectSimplex(points)
    plt.scatter(projected[:,0], projected[:,1], **kwargs)              
    # Leave some buffer around the triangle for vertex labels
    fig.gca().set_xlim(-0.2, 1.2)
    fig.gca().set_ylim(-0.2, 1.2)

    return fig    

def projectSimplex(points):
    """ 
    Project probabilities on the 3-simplex to a 2D triangle
    
    N points are given as N x 3 array
    """
    # Convert points one at a time
    tripts = np.zeros((points.shape[0],2))
    for idx in range(points.shape[0]):
        # Init to triangle centroid
        x = 1.0 / 2
        y = 1.0 / (2 * np.sqrt(3))
        # Vector 1 - bisect out of lower left vertex 
        p1 = points[idx, 0]
        x = x - (1.0 / np.sqrt(3)) * p1 * np.cos(np.pi / 6)
        y = y - (1.0 / np.sqrt(3)) * p1 * np.sin(np.pi / 6)
        # Vector 2 - bisect out of lower right vertex  
        p2 = points[idx, 1]  
        x = x + (1.0 / np.sqrt(3)) * p2 * np.cos(np.pi / 6)
        y = y - (1.0 / np.sqrt(3)) * p2 * np.sin(np.pi / 6)        
        # Vector 3 - bisect out of top vertex
        p3 = points[idx, 2]
        y = y + (1.0 / np.sqrt(3) * p3)
      
        tripts[idx,:] = (x,y)

    return tripts

# Define different colors for each label
cmap = CM.get_cmap('spectral')
norm = C.Normalize(vmin=0, vmax=len(labels))
c = range(len(labels))
# Do scatter plot
#fig = plotSimplex(testpoints, s=100, c=c, cmap=cmap, norm=norm)

probs = data[data.aa == "M"].set_index(["ff", "water", "aa"])
fig = plotSimplex(probs, s=100)
