#######################################################
# CALIBRATION
#######################################################
# FiveBars robot Calibration
import math
from math import sin,cos
import numpy
import robot

#~ from importlib import reload
#~ import calibration,robot
# RR kinematic functions
def f_5R(architecture,pose,command):
    [a11,a12,a21,a22,l1,l2,l3,l4] = architecture
    [x1,x2] = pose
    [q1,q2] = numpy.radians(command)
    f1 = (x1-a11-l1*numpy.cos(q1))**2 + (x2-a12-l1*numpy.sin(q1))**2 - l2**2
    f2 = (x1-a21-l4*numpy.cos(q2))**2 + (x2-a22-l4*numpy.sin(q2))**2 - l3**2
    return [f1,f2]

# Actuation of the robot in order to generate measures for calibration
def make_measurements(r5,commands,col='black',mar='*'):
    r5.actuate(commands[0])
    measures=[]
    print('   Taking measures ...')
    for q in commands:
        # print("q : ",q)
        r5.actuate(q)
        x = r5.measure_pose()
        r5.ax.plot([x[0]],[x[1]],color=col,marker=mar)
        # if syngularity met, don't add it to measures
        measures.append((x,q))
    r5.go_home()
    return measures

# Generation de commandes
def gen_command():
    commands = [[q1,q2] for q1 in range(0,45,5) for q2 in range(80,200,5)] + [[q1,q2] for q1 in range(5,50,5) for q2 in range(110,245,5) if (q1 != 45 and q2 != 240) ] + [[60,q2] for q2 in range(100,220,5)] + [[70,q2] for q2 in range(105,215,5)] + + [[q1,180] for q1 in range(-30,100,5)] +[[80,q2] for q2 in range(115,205,5)] + [[90,q2] for q2 in range(125,190,5)] + [[-30,q2] for q2 in range(110,180,5)] + [[-25,q2] for q2 in range(105,180,5)] + [[-20,q2] for q1 in range(-20,70,5) for q2 in range(100,180,5)] + [[70,q2] for q2 in range(105,180,5)] + [[80,q2] for q2 in range(115,180,5)] + [[85,q2] for q2 in range(120,180,5)] + [[90,q2] for q2 in range(125,180,5)] + [[95,q2] for q2 in range(130,180,5)]
    return commands

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


##################################################
# Execution
##################################################
# Definition of a 5R robot with approximate architecture
nominal_architecture = [-22.5,0,22.5,0,17.8,17.8,17.8,17.8]
r5 = robot.FiveBars(nominal_architecture,mode=0,seed=1)
#~ r5 = robot.TwoBars(nominal_architecture,mode=1,seed=1)

# Commands generation
commands = gen_command()
# Make measurements
measures = make_measurements(r5,commands)
# Calibration measure
calibrated_architecture = calibrate(f_5R,nominal_architecture,measures)
##################################################



##################################################
# Resultat de la calibration
# r5 = robot.FiveBars(nominal_architecture,mode=0,seed=1)
# nominal_architecture = [-22.5,0,22.5,0,17.8,17.8,17.8,17.8]
############### solution #############################
# >>> calibrated_architecture = calibration.calibrate(calibration.f_5R,nominal_architecture,measures)
#    Calibration processing ...
#    status :  Both `ftol` and `xtol` termination conditions are satisfied.
#    error :  115.55923428269708
#    result :  [-22.24157649   0.2888078   22.51328886  -0.15281802  17.58176891
#   17.80539455  17.51237713  17.59555565]
# architecture = [-22.24157649,0.2888078,22.51328886,-0.15281802,17.58176891,17.80539455,17.51237713,17.59555565]
##################################################
