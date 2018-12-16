"""Classes and functions to manipulate sets of boxes (pavings)

projected_centers_distance : box-distance function

Box : class representing a box

Paving : class representing a paving
"""
from __future__ import print_function

__authors__ = "A.Goldsztejn and C.Jermann"
__date__ = "01.09.2017"
__version__ = "2.1"
__copyright__ = "Copyright 2017-today, Universite de Nantes"

from matplotlib import pyplot
import math
import numpy
import rtree
from scipy import sparse
import struct

########################################################################
# Helper functions

def projected_centers_distance(box1,box2,dims):
    '''
    Function computing the distance between the projected centers of two boxes

    Parameters:
        box1,box2 - Two boxes in D dimensions
        dims - an array of dimension numbers (in 1..D) for the projection
    '''
    assert len(box1.vec) == len(box2.vec), "Cannot compute projected centers distance between boxes of different dimensions"
    assert 1<=min(dims) and max(dims)<=1+len(box1.vec)/2, "Invalid centers projection"
    dis = 0
    for d in dims:
        d -= 1
        dis += 0.25*((box1.vec[2*d]+box1.vec[2*d+1]) - (box2.vec[2*d]+box2.vec[2*d+1]))**2
    return math.sqrt(dis)

########################################################################
# Box class

class Box:
    '''
    Class Box holds the box coordinates and properties as computed by IBEX

    Members:
        idx - The box index ; should be unique
        vec - The box coordinates (x1min, x1max, x2min, x2max, ...)
        dat - optional data dictionary associated to the box (e.g. 'type' The box type (0=inner, 1=boundary, 2=unknown, 3=pending), 'parameters' The box certification parameters (None if type>1))

    Methods:
        __init__ - constructor
        __repr__ - printer
        __eq__ - equality comparator
        intersects - predicate
        draw2D - drawer
    '''

    def __init__(self,idx,vec,dat=None):
        '''
        Box constructor

        Parameters:
            idx - The box index ; should be unique
            vec - The box coordinates (x1min, x1max, x2min, x2max, ...)
            dat - optional data dictionary associated to the box (e.g. 'type' The box type (0=inner, 1=boundary, 2=unknown, 3=pending), 'parameters' The box certification parameters (None if type>1))
        '''
        self.idx = idx
        self.vec = vec
        self.dat = dat

    def __repr__(self):
        '''
        Box prettyprinter
        '''
        s = 'box ' + str(self.idx) + ' : ' + str(self.vec)
        for k in self.dat.keys():
            s += '  ' + str(k) + ':' + str(self.dat[k])
        return s

    def __eq__(self,other):
        '''
        Box equality comparison

        Parameter:
            other: the other box to compare to

        Return: True iff the boxes have the same bounds and associated data (though their indexes may differ)
        '''
        return self.vec==other.vec and self.dat==other.dat

    def intersects(self,bx):
        '''
        Box intersection predicate

        Parameter:
            bx - The box to check intersection with

        Return: True iff self and bx intersect
        '''
        assert len(self.vec)==len(bx.vec), 'Cannot compare boxes of different dimensions'
        if self.idx==bx.idx:
            return True
        for i in range(0,len(bx.vec)//2):
            [l1,u1] = self.vec[2*i:2*i+2]
            [l2,u2] = bx.vec[2*i:2*i+2]
            if not( (l1<=l2<=u1) or (l1<=u2<=u1) or (l2<=l1<=u2) or (l2<=u1<=u2) ):
                return False
        return True

    def draw2D(self,ax,x,y,ec='blue',fc='none'):
        '''
        Box 2D graphic drawer

        Parameters:
            ax - A matplotlib.axes.Axes (e.g., obtained from fig,ax=matplotlib.pyplot.subplots()) in which to display the box
            x - abscissa dimension (value: 1..n, with n=len(vec)//2)
            y - ordinate dimension (value: 1..n, with n=len(vec)//2)
            ec - color for the edges of the box
            fc - color for filling the box

        Return: nothing, but adds plot elements corresponding to the box to the passed artist
        '''
        assert 0<x<=len(self.vec)//2, "Incorrect dimension x"
        assert 0<y<=len(self.vec)//2, "Incorrect dimension y"
        x-=1 ; y-=1
        [xmin, xmax] = self.vec[2*x:2*x+2]
        [ymin, ymax] = self.vec[2*y:2*y+2]
        ax.add_artist(pyplot.Rectangle((xmin,ymin),xmax-xmin,ymax-ymin,edgecolor=ec,facecolor=fc))
        # alternative drawing method :
        # less efficient at construction, better at refreshing
        #~ ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color=ec)

########################################################################
# Paving class

class Paving:
    '''
    Class Paving holds a collection of boxes

    Members:
        boxes - the collection of boxes
        data - dictionary of associated data

    Methods:
        __init__ - constructor
        hull - bounding box
        boxes_intersecting - box selector
        subpaving - subpaving constructor
        load_mnf - IBEX MNF loader
        save_mnf - IBEX MNF saver
        draw2D - drawer
        adjacency_lists - grapher
        adjacency_matrix - grapher

    The paving also uses internally a RTree data structure to efficiently answer spatial indexing queries like "give all the boxes intersecting a given point/box" ; the corresponding member should not be used directly, but only through the provided methods.
    '''

    def __init__(self):
        '''
        Paving constructor, generate an empty paving
        '''
        self.boxes=[]
        self.data={}
        self.mnf_entete = "IBEX MANIFOLD FILE \0"
        self.mnf_version = 4

    def hull(self,d=None):
        '''
        Returns the hull, in dimensions d, of the paving content
        '''
        b = self._rtree.get_bounds()
        if d==None:
            res = b
        else:
            res = []
            for i in d:
                res.append(b[2*i-2])
                res.append(b[2*i-1])
        return res

    def boxes_intersecting(self,l,u=None,d=None):
        '''
        Boxes extractor

        Parameters:
            l - lowerbounds of the intersecting box
            u - upperbounds of the intersecting box (None if point)
            d - dimensions of the intersecting box (None if all dimensions)

        Return: the list of box indexes in the paving that intersect [l,u] in dimensions d
        '''
        if d==None:
            d=range(1,self.data['nvar']+1)
        if u==None:
            u=l
        b=self._rtree.get_bounds()
        for i in d:
            b[2*i-2],b[2*i-1] = l[i-1], u[i-1]
        return self._rtree.intersection(b)

    def subpaving(self,idx):
        '''
        Subpaving constructor

        Parameter:
            idx - Indexes of the boxes to include in the subpaving

        Return: a Paving composed of the boxes whose indexes are given in idx
        '''
        spv = Paving()
        # copying & updating global data
        for k in self.data.keys():
            spv.data[k] = self.data[k]
        spv.data['nin'] = 0
        spv.data['nbd'] = 0
        spv.data['nun'] = 0
        spv.data['npd'] = 0
        # copying boxes
        p = rtree.index.Property()
        p.dimension = self.data['nvar']
        spv._rtree = rtree.index.Index(properties=p,interleaved=False)
        for i in range(len(idx)):
            vec = self.boxes[idx[i]].vec
            typ = self.boxes[idx[i]].dat['type']
            # updating the box types counters
            if typ==0:    # inner
                spv.data['nin'] += 1
            elif typ==1:    # boundary
                spv.data['nbd'] += 1
            elif typ==2:     # unknown
                spv.data['nun'] += 1
            elif typ==3:     # pending
                spv.data['npd'] += 1
            par = self.boxes[idx[i]].dat['parameters']
            spv.boxes.append(Box(i,vec,{'type':typ,'parameters':par}))     # pushing the copied box in the subpaving
            spv._rtree.insert(i,vec)     # and in the rtree
        spv.data['nbox'] = spv.data['nin'] + spv.data['nbd'] + spv.data['nun'] + spv.data['npd']     # updating the total number of boxes
        return spv

    def load_mnf(self,filename,bxtyp=[0,1,2,3]):
        '''
        Paving loader from an IBEX manifold binary output file

        Parameter:
            filename - The name of the MNF file output by IBEX
            bxtyp - A list of box types of interest to be loaded from the MNF (default : [0,1,2,3] = all boxes)

        Return: nothing, but fills in the paving with boxes from the file
        '''
        with open(filename,'rb') as f:
            content = f.read()
            dic = {'pos':0}
            sizc = 1     # sizeof characters
            sizi = 4     # sizeof integers
            sizd = 8     # sizeof doubles

            def readstr(k=1):     # helper function to unpack strings
                pos = dic['pos']
                dic['pos'] += sizc*k
                return struct.unpack_from(str(k)+'s',content,pos)

            def readint(k=1):     # helper function to unpack (unsigned) integers
                pos = dic['pos']
                dic['pos'] += sizi*k
                return struct.unpack_from(str(k)+'I',content,pos)

            def readdbl(k=1):     # helper function to unpack (double precision) floats
                pos = dic['pos']
                dic['pos'] += sizd*k
                return struct.unpack_from(str(k)+'d',content,pos)

            # decoding preamble
            entete = readstr(20)[0].decode('utf-8')
            (version,) = readint()
            # checking version
            print(" Loading", str(entete), "version", version)
            assert (entete==self.mnf_entete and version<=self.mnf_version), 'bad version input file'
            # decoding global data
            (self.data['nvar'],) = readint()     # number of variables
            (self.data['nce'],) = readint()     # number of constraint equations
            (self.data['nci'],) = readint()     # number of constraint inequalities
            if version>=3: # MNF v3+ additional data extraction
                self.data['vars'] = [] # variable names
                for i in range(0,self.data['nvar']):
                    # reading char by char until a space is read
                    var = ""
                    ch=readstr(1)[0].decode('utf-8')
                    while ch!=" " and ch!="\0":
                        var = var + ch
                        ch=readstr(1)[0].decode('utf-8')
                    self.data['vars'].append(var)
            (self.data['sta'],) = readint()     # computation status
            (self.data['nin'],) = readint()     # number of inner boxes
            (self.data['nbd'],) = readint()     # number of boundary boxes
            (self.data['nun'],) = readint()     # number of unknown boxes
            (self.data['npd'],) = readint()     # number of pending boxes
            (self.data['time'],) = readdbl()     # computation time
            (self.data['ncel'],) = readint()     # number of cells
            # ~ print(self.data) ###################################################
            # decoding boxes
            self.data['nbox'] = self.data['nin'] + self.data['nbd'] + self.data['nun'] + self.data['npd']
            self.boxes=[]
            p = rtree.index.Property()
            p.dimension = self.data['nvar']
            self._rtree = rtree.index.Index(properties=p,interleaved=False)
            for i in range(self.data['nbox']):
                print("   box {}/{}\r".format(i+1,self.data['nbox']),end='\r')
                vec = readdbl(2*self.data['nvar'])     # box coordinates
                (typ,) = readint()    # box type
                if 0 < self.data['nce'] < self.data['nvar'] and ((version==1 and typ<2) or version>=2):
                    par = readint(self.data['nvar']-self.data['nce'])     # certificate parameters
                else:
                    par = None     # no certificate parameters
                if typ in bxtyp:
                    self.boxes.append(Box(i,vec,{'type':typ,'parameters':par}))     # pushing the decoded box in the paving
                    self._rtree.insert(i,vec)
            print("")
            # updating number of boxes depending on selected types
            if not 0 in bxtyp:
                self.data['nin']=0
            if not 1 in bxtyp:
                self.data['nbd']=0
            if not 2 in bxtyp:
                self.data['nun']=0
            if not 3 in bxtyp:
                self.data['npd']=0
            self.data['nbox'] = self.data['nin'] + self.data['nbd'] + self.data['nun'] + self.data['npd']

    def save_mnf(self,filename):
        '''
        Paving saver as an IBEX manifold binary file

        Parameter:
            filename - The name of the MNF file that can be input to IBEX

        Return: nothing, but fills in the file with boxes from the paving
        '''
        with open(filename,'wb') as f:
            content = bytearray()     # binary buffer

            def writestr(s):     # helper function to pack strings
                f.write(struct.pack(str(len(s))+'s',bytes(s,"utf-8")))

            def writeint(i):     # helper function to pack (unsigned) integers
                f.write(struct.pack('I',i))

            def writedbl(d):     # helper function to pack (double precision) floats
                f.write(struct.pack('d',d))

            # encoding preamble
            writestr(self.mnf_entete)
            writeint(2)
            print(" Writing", str(self.mnf_entete), "version", 2, " to ",filename)
            # encoding global data
            writeint(self.data['nvar'])     # number of variables
            writeint(self.data['nce'])     # number of constraint equations
            writeint(self.data['nci'])     # number of constraint inequalities
            writeint(self.data['sta'])     # computation status
            writeint(self.data['nin'])     # number of inner boxes
            writeint(self.data['nbd'])     # number of boundary boxes
            writeint(self.data['nun'])     # number of unknown boxes
            writeint(self.data['npd'])     # number of pending boxes
            writedbl(self.data['time'])     # computation time
            writeint(self.data['ncel'])     # number of cells
            # encoding boxes
            for b in self.boxes:
                for d in b.vec:     # box coordinates
                    writedbl(d)
                writeint(b.dat['type'])    # box type
                if 0 < self.data['nce'] < self.data['nvar'] and b.dat['type']<2:
                    for p in b.dat['parameters']:     # certificate parameters
                        writeint(p)

    def draw2D(self,ax,x,y,typ=False,ec='blue',fc='none'):
        '''
        Paving 2D graphic drawer

        Parameters:
            ax - A matplotlib.axes.Axes (e.g., obtained from fig,ax=matplotlib.pyplot.subplots()) in which to display the paving
            x - abscissa dimension (value: 1..n, with n=len(vec)//2)
            y - ordinate dimension (value: 1..n, with n=len(vec)//2)
            typ - draw box using colors representing their type (green=inner,blue=boundary,red=unknown,black=pending)
            ec - edge color
            fc - fill color

        Return: nothing, but adds plot elements corresponding to the box to the passed artist
        '''
        for bx in self.boxes:
            if typ==True:
                if bx.dat['type']==0:    # inner
                    ec='green'
                elif bx.dat['type']==1:    # boundary
                    ec='blue'
                elif bx.dat['type']==2:     # unknown
                    ec='red'
                elif bx.dat['type']==3:     # pending
                    ec='black'
                else:    # unrecognized type
                    ec='gray'
                bx.draw2D(ax,x,y,ec=ec)
            else:
                bx.draw2D(ax,x,y,ec=ec,fc=fc)

    def adjacency_lists(self):
        '''
        Paving grapher

        Return: a graph of overlapping boxes from the paving, in the form of a list of adjacency lists
        '''
        lst = []
        # generate intersections by rtree queries
        for bx in self.boxes:
            lst.append(list(self._rtree.intersection(bx.vec)))
        return lst

    def adjacency_matrix(self,weight=lambda x,y : 1):
        '''
        Paving grapher

        Parameter:
            weight - a function taking two boxes (the arc extremities) as parameters and returning a number to weight their connection

        Return: a graph of overlapping boxes from the paving, in the form of a scipy.sparse.csr_matrix
        '''
        assert(callable(weight))
        data = []
        row = []
        col = []
        # generate intersections by rtree queries
        for bx in self.boxes:
            l = list(self._rtree.intersection(bx.vec))     # query -> list of box indexes
            row.extend([bx.idx for k in range(len(l))])     # start node = box index
            col.extend(l)     # end nodes = query result
            data.extend([weight(bx,self.boxes[i]) for i in l])     # arc weights
        return sparse.csr_matrix((data, (row, col)), shape=(len(self.boxes), len(self.boxes)), dtype=numpy.float64)
