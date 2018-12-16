#######################################################
# TwoBars robot demo
import math
import numpy

#~ from importlib import reload
#~ import demo

##################################################
# Definition of a RR robot with approximate architecture
import robot

#~ nominal_architecture = [0,0,3,2]
#~ r2 = robot.TwoBars(nominal_architecture,seed=1,man_var=0.1,mes_var=0.01,eps_cmd=10)
#~ calibrated_architecture = r2.get_architecture() # fake calibration

###############################################################
# PATH PLANNING 

# define and display the obstacles
from matplotlib import pyplot

def put_obstacles(rob):
    obstacles = [ (-3, 2, 1) , (0, 3, 1) , (3, 2, 1) ]
    for o in obstacles:
        rob.ax.add_artist(pyplot.Circle((o[0],o[1]),o[2],color='.5',fill=False))
    rob.refresh()
    return obstacles

#~ obstacles = demo_planning.put_obstacles(r2)

# COMPUTATION OF R2 SINGULARITY&COLLISION-FREE CONFIGURATION SPACE USING IBEX

import paving

#~ import paving
#~ paving_RR = paving.Paving()
#~ paving_RR.load_mnf('RR-obs.mnf') # 2R with obstacles and without singularities
#~ neighborhood = paving_RR.adjacency_matrix()

from scipy.sparse.csgraph import dijkstra

# Path planning between two target boxes
def path_planner(pav,nei,Bori,Bdes):
    print('Building shortest path from ',Bori,' to ',Bdes)
    # construction of the shortest path
    path=[]
    dist,pred = dijkstra(nei,return_predecessors=True,indices=Bori)
    # checking existence of a compatible destination box
    if dist[Bdes]!=numpy.inf:
        print('   shortest path found!')
        # collecting the path from the predecessors tree
        path=[Bdes]
        while path[0]!=Bori:
            path.insert(0,pred[path[0]])
    else:
        print('   impossible to connect boxes!')
    return path

Xori,Xdes = [-3,0.5], [3,0.5]
r2.ax.plot([Xori[0],Xdes[0]],[Xori[1],Xdes[1]],linestyle=':',marker='*',color='.3')
r2.refresh()
ori = list(paving_RR.boxes_intersecting(Xori,d=[1,2]))
des = list(paving_RR.boxes_intersecting(Xdes,d=[1,2]))
shortest_path=demo_planning.path_planner(paving_RR,neighborhood,ori[0],des[0])


#~ from scipy.sparse.csgraph import connected_components
#~ nc,l = connected_components(neighborhood, directed=False)
#~ for b in [i for i in range(len(paving_RR.boxes)) if l[i]==l[ori[0]]]:
    #~ paving_RR.boxes[b].draw2D(r2.ax,1,2)
#~ for b in [i for i in range(len(paving_RR.boxes)) if l[i]!=l[ori[0]]]:
    #~ paving_RR.boxes[b].draw2D(r2.ax,1,2,ec='magenta')
#~ r2.refresh()
#~ del r2
#~ r2 = robot.TwoBars(nominal_architecture,seed=1,man_var=0.1,mes_var=0.01,eps_cmd=10)
#~ obstacles = demo_planning.put_obstacles(r2)


Xori,Xdes = [-4.5,1.5], [4.5,1.5]
r2.ax.plot([Xori[0],Xdes[0]],[Xori[1],Xdes[1]],linestyle=':',marker='*',markersize=10,color='.3')
r2.refresh()
ori = list(paving_RR.boxes_intersecting(Xori,d=[1,2]))
des = list(paving_RR.boxes_intersecting(Xdes,d=[1,2]))
#~ shortest_path=demo_planning.path_planner(paving_RR,neighborhood,ori[0],des[0])


# Display a box path and actuate the robot along it
def display_path(rob,pav,spath,bcol='yellow',rcol='cyan'):
    def midq(b):
        return numpy.degrees([(b.vec[4]+b.vec[5])/2,(b.vec[6]+b.vec[7])/2])
    if spath!=[]:
        rob.pen_up()
        rob.actuate(midq(pav.boxes[spath[0]]))
        rob.pen_down(rcol)
        for b in spath:
            if bcol!=None:
                pav.boxes[b].draw2D(rob.ax,1,2,ec=bcol,fc=bcol)
            rob.actuate(midq(pav.boxes[b]))
        rob.pen_up()
        rob.go_home()


#~ demo_planning.display_path(r2,paving_RR,shortest_path)

#~ neighborhood = paving_RR.adjacency_matrix(weight=lambda x,y : paving.projected_centers_distance(x,y,[1,2]))
