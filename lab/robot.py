'''Classes and functions to manipulate virtual robots

Robot : Base abstract class representing virtual robots

FiveBars : Class representing a five-bars (R_RRRR_) robot (derived from Robot)

TwoBars : Class representing a two-bars (R_R_) robot (derived from Robot)
'''
from __future__ import print_function

__authors__ = "A.Goldsztejn and C.Jermann"
__date__ = "01.09.2017"
__version__ = "1.0"
__copyright__ = "Copyright 2017-today, Universite de Nantes"

import random
import numpy
from matplotlib import pyplot,patches

########################################################################
# Helper functions

def linear_discretization(a, b, eps):
    '''Computation of linear steps of norm <= eps between a and b.

        Arguments:
        a -- a point (numpy.array)
        b -- a vector/delta (numpy.array)
        eps -- the maximal norm

        Return:
        A list of equally distributed aligned points between a (excluded) and b (included)
    '''
    points=[]
    a=numpy.array(a)
    b=numpy.array(b)
    if max(abs(a-b))!=0:
        u=(b-a)/numpy.ceil(numpy.linalg.norm(b-a)/eps)
        count=numpy.linalg.norm(b-a)/numpy.linalg.norm(u)
        for i in range(int(count)-1):
            points.append(a + (i+1)*u)
        points.append(b)
    return points

def circles_intersections(c01,c02,r0,c11,c12,r1):
    '''Computation of the intersection of two circles

        The circle equations are:
        1) (x-c01)^2+(y-c02)^2=r0^2
        2) (x-c11)^2+(y-c12)^2=r1^

        Arguments:
        c01, ... - circles coordinates (see above)

        Return:
        A numpy.array of the coordinates of the intersections

        Remark: returns an empty numpy.array if the two circles do not intersect or completely overlap
        '''
    d = numpy.sqrt((c11 - c01)**2 + (c12 - c02)**2)
    #
    rt0 = r0/d
    rt1 = r1/d
    ct11 = (c11 - c01)/d
    ct12 = (c12 - c02)/d
    D = -(-1 + rt0 - rt1) * (1 + rt0 - rt1) * (-1 + rt0 + rt1) * (1 + rt0 + rt1)
    #
    if D<0:
        return numpy.array([])
    else:
        sD=numpy.sqrt(D)
    #
    x1 = c01 + d/2 * ( ct12 * sD + ct11 * (1 + rt0**2 - rt1**2))
    y1 = c02 + d/2 * (-ct11 * sD + ct12 * (1 + rt0**2 - rt1**2))
    #
    x2 = c01 + d/2 * (-ct12 * sD + ct11 * (1 + rt0**2 - rt1**2))
    y2 = c02 + d/2 * ( ct11 * sD + ct12 * (1 + rt0**2 - rt1**2))
    #
    if D>0:
        return numpy.array([[x1,y1],[x2,y2]])
    else:
        return numpy.array([[x1,y1]])

########################################################################
# Virtual robot class (generic)

class Robot:
    '''A controlable virtual robot.

    Members:
        ax - The matplotlib.axes.Axes in which is displayed the robot ; can be used to overlay pavings and other figures

    Methods:
        __init__ - Constructor
        __del__ - Destructor
        pen_down - Tracing activation
        pen_up - Tracing deactivation
        actuate - Actuation with absolute command
        actuate_rel - Actuation with relative command
        go_home - Actuation to home command
        measure_pose - Measurement of the current pose
        measure_command - Measurement of the current command
        refresh - Refreshes the robot display

    Remarks:
    * This class is just an interface, it must be specialized for actual robots which must at least overload the kinematic models (see FiveBars below).
    * The public member can be accessed and changed anytime without harming the robot behavior, allowing the superimposition of additional graphic elements
    '''

    def __init__(self, architecture, mode, seed, man_var, mes_var, home_cmd, eps_cmd):
        '''Robot constructor

        Parameters:
            architecture - A list of the robot architecture parameters (attach points, bar lengths, etc.)
            mode - The initial assembly mode
            seed - The tolerances seed
            man_var - The manufacturing tolerance variance
            mes_var - The measurement captors variance
            home_cmd - The initial robot command
            eps_cmd - The actuation maximal step size

        A robot may have manufacturing imprecisions, yielding perturbated architecture parameters and perturbated measurement captors. These are modeled using a random Gaussian distribution around the nominal/actual values parameterized by the provided seed, the manufacturing variance man_var, and the measurement variance mes_var.
        '''
        random.seed(a=seed)
        #~ self._wmode=[ 1+random.random() for a in architecture ]
        self._architecture = [ random.gauss(p, man_var) for p in architecture]  # perturbated architecture generation
        #~ self._architecture = [ (k+random.gauss(a, man_var))/k for (a,k) in zip(architecture,self._wmode) ]    # encrypted perturbated architecture generation
        self._mes_var = mes_var
        self._home_cmd=home_cmd
        self._eps_cmd=eps_cmd
        self._cmd = numpy.array(self._home_cmd)    # conversion to numpy.array
        self._mode = mode
        self._pos = self._direct_kinematic_model(self._cmd,self._mode)
        self._pen = False
        self._fig, self.ax = pyplot.subplots()
        self.ax.set_aspect('equal')
        self._draw_backup=self._draw_robot()
        self._draw_workspace()
        self._fig.show()

    def __del__(self):
        '''Robot destructor'''
        pyplot.close(self._fig)

    def pen_down(self,color='green',width=1):
        '''Put the pen down so future trajectories will be traced.'''
        self._pen = True
        self._pen_color=color
        self._pen_width=width
        self._redraw_robot()
        self._fig.canvas.draw()

    def pen_up(self):
        '''Put the pen up so future trajectories will not be traced.'''
        self._pen = False
        self._redraw_robot()
        self._fig.canvas.draw()

    def actuate(self, cmd):
        '''Actuate the robot to change its current commands to cmd.'''
        discretized_cmds=linear_discretization(self._cmd,cmd,self._eps_cmd)
        for q in discretized_cmds:
            try:
                x=self._direct_kinematic_model(q,self._mode)
                if self._pen:
                    self.ax.plot([self._pos[0],x[0]],[self._pos[1],x[1]],color=self._pen_color,linewidth=self._pen_width)
                self._pos,self._cmd=x,q
                self._redraw_robot()
                self._fig.canvas.draw()
            except ValueError as e:
                print('!!! ERROR: ',e)
                break
        self._redraw_robot()
        self._fig.canvas.draw()

    def actuate_rel(self, cmd_delta):
        '''Actuate the robot to change its current commands by cmd_delta.'''
        cmd_delta = numpy.array(cmd_delta)    # conversion to numpy.array
        self.actuate(self._cmd+cmd_delta)

    def go_home(self):
        '''Actuate the robot so it returns to its home command'''
        self.actuate(self._home_cmd)

    def measure_command(self):
        '''Return a measurement of the current command of the robot'''
        return self._cmd

    def measure_pose(self):
        '''Return a measurement of the current pose of the robot. since the captors may be imprecise, this measurement may be noisy (see Robot.__init__)'''
        return numpy.array([ random.gauss(p, self._mes_var) for p in self._pos])

    def refresh(self):
        '''Refresh robot display'''
        self._redraw_robot()
        self._fig.canvas.draw()

    def _direct_kinematic_model(self, cmd, mode):
        # returns the poses corresopnding to the given command and assembly mode
        # raise 'singularity' error is the number of solutions is not maximal
        raise NotImplementedError('Call to virtual _direct_kinematic_model: this method should be overloaded.')

    def _redraw_robot(self):
        # refreshes the current graphic by removing the old pose and plotting the new one
        for obj in self._draw_backup:
            obj.remove()
        self._draw_backup=self._draw_robot()

    def _draw_robot(self):
        # draws to robot current pose
        raise NotImplementedError('Call to virtual _draw_robot: this method should be overloaded.')

    def _draw_workspace(self):
        raise NotImplementedError('Call to virtual _draw_workspace: this method should be overleaded.')

########################################################################
# Fivebars robot class (specific)
# Dextar : [-22.5,0,22.5,0,17.8,17.8,17.8,17.8]

class FiveBars(Robot):
    '''The virtual 5-bars parallel robot

    Its _architecture is defined using 8 parameters [a11,a12,a21,a22,a31,a32,a41,a42] representing:
    a11,a12 -- coordinates of first anchor point
    a21,a22 -- coordinates of second anchor point
    a31,a32 -- lengths of arm and forearm between first anchor point and effector
    a41,a42 -- lengths of arm and forearm between second anchor point and effector

    It accepts two assembly modes (values : 0, 1)

    Its commands are the two angles [q1,q2] (in degrees) formed by each arm with the horizon

    Its pose is defined by the two Cartesian coordinates [x1,x2] of its end-effector
    '''

    def __init__(self, architecture, mode, seed, man_var=0.2, mes_var=0.01, home_cmd=[0,180], eps_cmd=2):
        '''FiveBars constructor (see Robot.__init__)'''
        Robot.__init__(self,architecture,mode,seed,man_var,mes_var,home_cmd,eps_cmd)

    def _direct_kinematic_model(self, cmd, mode):
        #~ [x1,y1,x2,y2,L1,l1,L2,l2] = [ a*k-k for (a,k) in zip(self._architecture,self._wmode) ]
        [x1,y1,x2,y2,L1,l1,L2,l2] = self._architecture
        [q1,q2] = numpy.radians(cmd)
        # elbow1 pose
        e1x, e1y = x1+L1*numpy.cos(q1), y1+L1*numpy.sin(q1)
        # elbow2 pose
        e2x, e2y = x2+L2*numpy.cos(q2), y2+L2*numpy.sin(q2)
        # effector poses
        poses = circles_intersections(e1x,e1y,l1,e2x,e2y,l2)
        if len(poses)!=2:
            raise ValueError('singularity met')
        return poses[mode]

    def _draw_robot(self):
        #~ [x1,y1,x2,y2,L1,l1,L2,l2] = [ a*k-k for (a,k) in zip(self._architecture,self._wmode) ]
        [x1,y1,x2,y2,L1,l1,L2,l2] = self._architecture
        # robot state
        [q1,q2] = numpy.radians(self._cmd)
        [X,Y] = self._pos
        # elbow1 pose
        e1x, e1y = x1+L1*numpy.cos(q1), y1+L1*numpy.sin(q1)
        # elbow2 pose
        e2x, e2y = x2+L2*numpy.cos(q2), y2+L2*numpy.sin(q2)
        # plotting the robot pose
        draw1=self.ax.plot([x1,e1x,X,e2x,x2],[y1,e1y,Y,e2y,y2],color='black',linestyle='solid',label='',marker='o')
        if self._pen:
            m,w='v',self._pen_color
        else:
            m,w='^','w'
        draw2=self.ax.plot([X],[Y],color=w,linestyle='solid',marker=m,markersize=9)
        return draw1+draw2

    def _draw_workspace(self):
        #~ [x1,y1,x2,y2,L1,l1,L2,l2] = [ a*k-k for (a,k) in zip(self._architecture,self._wmode) ]
        [x1,y1,x2,y2,L1,l1,L2,l2] = self._architecture
        xmin, xmax  = min(x1,max(x1-L1-l1,x2-L2-l2)), max(x2,min(x1+L1+l1,x2+L2+l2))
        xran = xmax-xmin
        ymin, ymax  = max(y1-L1-l1,y2-L2-l2), min(y1+L1+l1,y2+L2+l2)
        yran = ymax-ymin
        self.ax.axis([xmin-xran/10, xmax+xran/10, ymin-yran/10, ymax+yran/10]) # TODO : make a better range to focus on the actual workspace
        C1 = pyplot.Circle((x1,y1),L1+l1,color='.5',fill=False)
        c1 = pyplot.Circle((x1,y1),numpy.abs(L1-l1),color='.5',fill=False)
        C2 = pyplot.Circle((x2,y2),L2+l2,color='.5',fill=False)
        c2 = pyplot.Circle((x2,y2),numpy.abs(L2-l2),color='.5',fill=False)
        self.ax.add_artist(C1)
        self.ax.add_artist(c1)
        self.ax.add_artist(C2)
        self.ax.add_artist(c2)

########################################################################
# Twobars robot class (specific)
# RR : [0,0,20,15]

class TwoBars(Robot):
    '''The virtual 2-bars serial robot

    Its _architecture is defined using 4 parameters [a11,a12,a21,a22] representing:
    a11,a12 -- coordinates of the anchor point
    a21,a22 -- lengths of arm and forearm between anchor point and effector

    It accepts no assembly modes

    Its commands are the two angles [q1,q2] (in degrees) formed by the arm (q1) and fore-arm (q2) with the horizon

    Its pose is defined by the two Cartesian coordinates [x1,x2] of its end-effector ; X2 must be positive
    '''

    def __init__(self, architecture, seed=1, man_var=1.0, mes_var=0.1, home_cmd=[90,0], eps_cmd=2):
        '''TwoBars constructor (see Robot.__init__)'''
        Robot.__init__(self,architecture,0,seed,man_var,mes_var,home_cmd,eps_cmd)

    def get_architecture(self):
        return self._architecture

    def _direct_kinematic_model(self, cmd, mode):
        #~ [x,y,L,l] = [ a*k-k for (a,k) in zip(self._architecture,self._wmode) ]
        [x,y,L,l] = self._architecture
        [q1,q2] = numpy.radians(cmd)
        # elbow pose
        ex, ey = x+L*numpy.cos(q1), y+L*numpy.sin(q1)
        # effector pose
        pose = [ex+l*numpy.cos(q2), ey+l*numpy.sin(q2)]
        if pose[1]<y:
            raise ValueError('out of workspace')
        return pose

    def _draw_robot(self):
        #~ [x,y,L,l] = [ a*k-k for (a,k) in zip(self._architecture,self._wmode) ]
        [x,y,L,l] = self._architecture
        # robot state
        [q1,q2] = numpy.radians(self._cmd)
        [X,Y] = self._pos
        # elbow pose
        ex, ey = x+L*numpy.cos(q1), y+L*numpy.sin(q1)
        # plotting the robot pose
        draw1=self.ax.plot([x,ex,X],[y,ey,Y],color='black',linestyle='solid',label='',marker='o')
        if self._pen:
            m,c='v',self._pen_color
        else:
            m,c='^','w'
        draw2=self.ax.plot([X],[Y],color=c,marker=m,markersize=9)
        return draw1+draw2

    def _draw_workspace(self):
        #~ [x,y,L,l] = [ a*k-k for (a,k) in zip(self._architecture,self._wmode) ]
        [x,y,L,l] = self._architecture
        xmin, xmax  = x-L-l, x+L+l
        xran = xmax-xmin
        ymin, ymax  = y-min(L,l), y+L+l
        yran = ymax-ymin
        self.ax.axis([xmin-xran/10, xmax+xran/10, ymin-yran/10, ymax+yran/10])
        C1 = patches.Arc((x,y),2*(L+l),2*(L+l),0,0,180,color='.5',fill=False,linestyle=':',linewidth=0.5)
        c1 = patches.Arc((x,y),2*numpy.abs(L-l),2*numpy.abs(L-l),0,0,180,color='.5',fill=False,linestyle=':',linewidth=0.5)
        self.ax.add_artist(C1)
        self.ax.add_artist(c1)
        self.ax.plot([x-L-l,x-numpy.abs(L-l)],[y,y],color='.5',linestyle=':',linewidth=0.5)
        self.ax.plot([x+numpy.abs(L-l),x+L+l],[y,y],color='.5',linestyle=':',linewidth=0.5)
