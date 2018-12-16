#######################################################
# TwoBars robot demo
import math
import numpy
import robot
# import calibration


##################################################
# Definition of a RR robot with approximate architecture
# robot creation routine
def create_5R():
    architecture = [-22.24157649,0.2888078,22.51328886,-0.15281802,17.58176891,17.80539455,17.51237713,17.59555565]
    r5 = robot.FiveBars(architecture,mode=0,seed=1)
    # r5 = robot.FiveBars(architecture,mode=1,seed=1)
    return r5

###############################################################
# PATH FOLLOWING
# Hypothesis : the calibrated architecture has been computed and is used for all computations in this part

# Definition of the target trajectory
def find_Pose_time(t):
    return (0+5*numpy.cos(t),-20+5*numpy.sin(t))

# Construction of discretized target path
def discretize(r5,trajectory,tmin,tmax,steps):
    target_path = [trajectory(t) for t in numpy.linspace(tmin,tmax,steps)]
    r5.ax.plot([x[0] for x in target_path],[x[1] for x in target_path],color='blue',linestyle=':',marker='+')
    r5.refresh()
    return target_path

# RR kinematic functions
def f_5R(architecture,pose,command):
    [a11,a12,a21,a22,l1,l2,l3,l4] = architecture
    [x1,x2] = pose
    [q1,q2] = numpy.radians(command)
    f1 = (x1-a11-l1*numpy.cos(q1))**2 + (x2-a12-l1*numpy.sin(q1))**2 - l2**2
    f2 = (x1-a21-l4*numpy.cos(q2))**2 + (x2-a22-l4*numpy.sin(q2))**2 - l3**2
    return [f1,f2]

from scipy.optimize import root

# continuation procedure
def continuation(function_xq,target_path,q0):
    qs=[]
    qk=q0
    for xk in target_path:
        res = root(lambda q: function_xq(xk,q), qk)
        if not res.success:
            print(res.message)
            break
        qk = res.x
        qs.append(qk)
    return qs

def draw_path(r5,commands,col='blue'):
    r5.pen_up()
    r5.go_home()
    r5.actuate(commands[0])
    r5.pen_down(col)
    real_path = []
    for q in commands:
        r5.actuate(q)
        real_path.append(numpy.array(r5.measure_pose()))
    r5.pen_up()
    r5.go_home()
    return real_path

# distance function used to measure error
def dist(x,y):
    return numpy.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2)

##################################################
# Execution
##################################################
architecture = [-22.24157649,0.2888078,22.51328886,-0.15281802,17.58176891,17.80539455,17.51237713,17.59555565]
r5 = create_5R()
target_path = discretize(r5,find_Pose_time,0,2*numpy.pi,50)
x0 = target_path[0]
# Resolution avec IBEX pour identifier qo
 # IBEX solution 1 # sol n°0 = ([-0.3540031261646856, -0.3540031261646834] ; [4.704071811846806, 4.704071811846809]) [inner]
 # IBEX solution 2 # sol n°1 = ([-0.3540031261646856, -0.3540031261646834] ; [3.27468684929581, 3.274686849295813]) [inner]
 # IBEX solution 3 # sol n°2 = ([-0.9262925380081022, -0.9262925380081004] ; [4.704071811846806, 4.704071811846809]) [inner]
 # IBEX solution 4 # sol n°3 = ([-0.9262925380081022, -0.9262925380081004] ; [3.27468684929581, 3.274686849295813]) [inner]
q0 = numpy.degrees([-0.3540031261646856,4.704071811846806])
commands = continuation(lambda x,q:f_5R(architecture,x,q),target_path,q0)
real_path = draw_path(r5,commands,col='red')

# q0 = numpy.degrees([-0.3540031261646856,3.27468684929581])
# commands = continuation(lambda x,q:f_5R(architecture,x,q),target_path,q0)
# real_path = draw_path(r5,commands,col='red')

# Cas avec singularite
# q0 = numpy.degrees([-0.9262925380081022,4.704071811846806])
# commands = continuation(lambda x,q:f_5R(architecture,x,q),target_path,q0)
# real_path = draw_path(r5,commands,col='red')
#
# q0 = numpy.degrees([-0.9262925380081022,3.27468684929581])
# commands = continuation(lambda x,q:f_5R(architecture,x,q),target_path,q0)
# real_path = draw_path(r5,commands,col='red')
##################################################

#max([demo.dist(x,rx) for (x,rx) in zip(target_path,real_path)])

#~ commands = demo.continuation(lambda x,q:demo.f_RR(calibrated_architecture,x,q),target_path,q0)
#~ real_path = demo.draw_path(r5,commands)
#~ max([demo.dist(x,rx) for (x,rx) in zip(target_path,real_path)])
