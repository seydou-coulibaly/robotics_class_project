#######################################################
# TwoBars robot demo
import math
import numpy

#~ from importlib import reload
#~ import demo

##################################################
# Definition of a RR robot with approximate architecture
import robot

# robot creation routine
def create_2R(nominal_architecture,manufacture_tolerance=0,measurement_precision=0):
    r2 = robot.TwoBars(nominal_architecture,seed=1,man_var=manufacture_tolerance/2,mes_var=measurement_precision/2,eps_cmd=10)
    return r2

#~ nominal_architecture = [0,0,3,2]
#~ r2 = demo.create_2R(nominal_architecture)
#~ r2 = demo.create_2R(nominal_architecture,0.2,0.02)


##################################################
# CALIBRATION

# RR kinematic functions
def f_RR(architecture,pose,command):
    [a1,a2,a3,a4] = architecture
    [x1,x2] = pose
    [q1,q2] = numpy.radians(command)
    f1 = a1+a3*numpy.cos(q1)+a4*numpy.cos(q1+q2) - x1
    f2 = a2+a3*numpy.sin(q1)+a4*numpy.sin(q1+q2) - x2
    return [f1,f2]

# Actuation of the robot in order to generate measures for calibration
def make_measurements(r2,commands,col='black',mar='*'):
    r2.actuate(commands[0])
    #~ r2.pen_down()
    measures=[]
    print('   Taking measures ...')
    for q in commands:
        r2.actuate(q)
        x = r2.measure_pose()
        r2.ax.plot([x[0]],[x[1]],color=col,marker=mar)
        measures.append((x,q))
    #~ r2.pen_up()
    r2.go_home()
    return measures

#~ commands = [[q,q] for q in range(0,1,1)]
#~ commands = [[q1,q2] for q1 in range(0,181,45) for q2 in range(-q1,181-q1,30)]
#~ measures = demo.make_measurements(r2,commands)

# calibration from measurements
from scipy.optimize import least_squares

def calibrate(kinematic_functions,nominal_architecture,measures):
    # error function
    def errors(a):
        err=[]
        for (x,q) in measures:
            for fi in kinematic_functions(a,x,q):
                err.append(fi)
        return err
    print('   Calibration processing ...')
    sol = least_squares(errors,nominal_architecture)
    print('   status : ',sol.message)
    print('   error : ',sol.cost)
    print('   result : ',sol.x)
    return sol.x

#~ calibrated_architecture = demo.calibrate(demo.f_RR,nominal_architecture,measures)






###############################################################
# PATH FOLLOWING 
# Hypothesis : the calibrated architecture has been computed and is used for all computations in this part

# Definition of the target trajectory
def lemniscate(t):
    cos2 = 1+numpy.cos(t)**2
    d1 = 2/cos2
    d2 = 2/cos2
    return (2+d1*numpy.sin(t),2+d2*numpy.cos(t)*numpy.sin(t))
    
# Construction of discretized target path
def discretize(r2,trajectory,tmin,tmax,steps):
    target_path = [trajectory(t) for t in numpy.linspace(tmin,tmax,steps)]
    r2.ax.plot([x[0] for x in target_path],[x[1] for x in target_path],color='blue',linestyle=':',marker='+')
    r2.refresh()
    return target_path

#~ import numpy
#~ target_path = demo.discretize(r2,demo.lemniscate,0,2*numpy.pi,50)
#~ x0 = target_path[0]
#~ q0 = numpy.degrees([0.097389470147109,2.095906191442448]) # IBEX solution 1
#~ q0 = numpy.degrees([1.464922832974473,-2.095906191442448]) # IBEX solution 2
#~ r2.actuate(q0)
    
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
        

#~ commands = demo.continuation(lambda x,q:demo.f_RR(nominal_architecture,x,q),target_path,q0)

def draw_path(r2,commands,col='blue'):
    r2.pen_up()
    r2.go_home()
    r2.actuate(commands[0])
    r2.pen_down(col)
    real_path = []
    for q in commands:
        r2.actuate(q)
        real_path.append(numpy.array(r2.measure_pose()))
    r2.pen_up()
    r2.go_home()
    return real_path

#~ real_path = demo.draw_path(r2,commands,col='red')

# distance function used to measure error
def dist(x,y):
    return numpy.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2)

#~ max([demo.dist(x,rx) for (x,rx) in zip(target_path,real_path)])

#~ commands = demo.continuation(lambda x,q:demo.f_RR(calibrated_architecture,x,q),target_path,q0)
#~ real_path = demo.draw_path(r2,commands)
#~ max([demo.dist(x,rx) for (x,rx) in zip(target_path,real_path)])
