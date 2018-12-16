import paving
p = paving.Paving()
p.load_mnf("simple.mnf")
from matplotlib import pyplot
fig1,ax1 = pyplot.subplots()
p.draw2D(ax1,1,2)
ax1.axis(p.hull([1,2]))
fig1.show()
# m = p.adjacency_matrix()
# from scipy.sparse.csgraph import connected_components
# nc,l = connected_components(m, directed=False)
# sp0 = p.subpaving([i for i in range(len(p.boxes)) if l[i]==0])
# sp0.draw2D(ax1,1,2,ec=None,fc='yellow')
#
# import robot
# seed = 0
# r=robot.FiveBars([-22.5, 0, 22.5, 0, 17.8, 17.8, 17.8,17.8],0,seed)
# r.measure_pose()
# r.pen_up()
# r.go_home()
# r.pen_down()
