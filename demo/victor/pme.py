import copy
import itertools
import os
import re

def minus(op):
    if op==" - ":
        return " + "
    else: return " - "
def stradd(a,b):
    at = a.split("_")[0]; bt = b.split("_")[0]
    op = " + "
    if a[0]=="-":
        a = a[1:]; op = minus(op)
    if b[0]=="-":
        b = b[1:]; op = minus(op)
    if at=="0": return b
    elif bt=="0": return a
    else: return a+op+b
def strsub(a,b):
    at = a.split("_")[0]; bt = b.split("_")[0]
    if at=="0": return b
    elif bt=="0": return a
    else: return a+" - "+b
def strmul(a,b):
    at = a.split("_")[0]; bt = b.split("_")[0]
    if at=="0" or bt=="0": return "0"
    elif at=="1" or at=="I": return b
    elif bt=="1" or bt=="I": return a
    elif at=="-1" or at=="-I": return "-"+b
    elif bt=="-1" or bt=="-I": return "-"+a
    else: return a+" * "+b

class Smatrix():
    def __init__(self,isize,jsize):
        self.isize = isize; self.jsize = jsize
        self.array = [ [ "0" for j in range(jsize) ] for i in range(isize) ]
        if isize>1:
            self.vletters = [ "m" for i in range(self.isize) ]
            self.vletters[0] = "t"; self.vletters[self.isize-1] = "b"
        else:
            self.vletters = [ "" for i in range(self.isize) ]
        if jsize>1:
            self.hletters = [ "m" for j in range(self.jsize) ]
            self.hletters[0] = "l"; self.hletters[self.jsize-1] = "r"
        else:
            self.hletters = [ "" for j in range(self.jsize) ]
    def __getitem__(self,idx):
        if idx[0]<0 or idx[0]>=len(self.array):
            print "ERROR i index %d out of range %d" % \
                  (idx[0],len(self.array)-1)
        if idx[1]<0 or idx[1]>=len(self.array[idx[0]]):
            print "ERROR j index %d out of range %d" % \
                  (idx[1],len(self.array[idx[0]])-1)
        return self.array[idx[0]][idx[1]]
    def __setitem__(self,idx,val):
        self.array[idx[0]][idx[1]] = val
    def multisub(self,iset):
        # iset is an array of subset expressions: [ [0,0], [1,0] ]
        s = []
        for i in iset:
            s.append(self[(i)])
        return s
    def setscalar(self,value):
        for k in range(self.isize):
            self[(k,k)] = value
        self.annotated()
    def isscalar(self):
        r = True
        for i in range(self.isize):
            for j in range(self.jsize):
                if i!=j:
                    r = r and self[(i,j)]!="0"
        return r
    def setfull(self,value):
        for i in range(self.isize):
            for j in range(self.jsize):
                self[(i,j)] = value
        self.annotate()
    def setupper(self,value):
        for i in range(self.isize):
            for j in range(self.jsize):
                if j>=i: self[(i,j)] = value
        self.annotate()
    def setunitupper(self,value):
        for i in range(self.isize):
            for j in range(self.jsize):
                if j==i and i>0 and i<self.isize-1 :
                    self[(i,j)] = "1"
                elif j>=i:
                    self[(i,j)] = value
        self.annotate()
    def setstrictupper(self,value):
        for i in range(self.isize):
            for j in range(self.jsize):
                if j>i or (j==i and (i==0 or i==self.isize-1) ):
                    self[(i,j)] = value
        self.annotate()
    def setlower(self,value):
        for i in range(self.isize):
            for j in range(self.jsize):
                if j<=i: self[(i,j)] = value
        self.annotate()
    def setlowerbi(self,value):
        for i in range(self.isize):
            for j in range(self.jsize):
                if j==i or j==i-1: self[(i,j)] = value
        self.annotate()
    def setI(self):
        for i in range(self.isize):
            for j in range(self.jsize):
                if j==i:
                    self[(i,j)] = "I"
        self.annotate()
    def setJ(self):
        for i in range(self.isize):
            for j in range(self.jsize):
                if j==i:
                    if i==0 or i==self.isize-1: self[(i,j)] = "J"
                elif i==j+1: self[(i,j)] = "j"
        self.annotate()
    def annotate(self):
        for i in range(self.isize):
            for j in range(self.jsize):
                e = self[(i,j)]; x = self.vletters[i]+self.hletters[j]
                if e!="0" and e!="1" and x!="":
                    self[(i,j)] = e+"_"+x
    def annotated(self):
        for i in range(self.isize):
            for j in range(self.jsize):
                e = self[(i,j)]; x = self.hletters[j]
                if e!="0" and e!="I" and x!="":
                    self[(i,j)] = e+"_"+x
    def __repr__(self):
        return repr(self.array)
    def __str__(self):
        if self.isize==1:
            s = "[ "
            for c in self.array[0]:
                if c=="0":
                    s += "_,"
                else: s += c+"," 
            s += " ]"
        else:
            s = "[ "
            for r in self.array:
                s += "[ "
                for c in r:
                    if c=="0":
                        s += "_,"
                    else: s += c+"," 
                s += " ], "
            s += "]"
        return s
    def t(self):
        r = Smatrix(self.jsize,self.isize)
        for i in range(self.isize):
            for j in range(self.jsize):
                r[(j,i)] = "T("+self[(i,j)]+")"
        return r
    def __add__(self,other):
        if self.isize!=other.isize or self.jsize!=other.jsize:
            print "Incompatible matrices"; return 0
        r = Smatrix(self.isize,self.jsize)
        for i in range(self.isize):
            for j in range(self.jsize):
                r[(i,j)] = stradd( self[(i,j)], other[(i,j)] )
        return r
    def __sub__(self,other):
        if self.isize!=other.isize or self.jsize!=other.jsize:
            print "Incompatible matrices"; return 0
        r = Smatrix(self.isize,self.jsize)
        for i in range(self.isize):
            for j in range(self.jsize):
                r[(i,j)] = strsub( self[(i,j)], other[(i,j)] )
        return r
    def __mul__(self,other):
        M = self.isize; N = other.jsize; K1 = self.jsize; K2 = other.isize
        if K1!=K2:
            print "Incompatible matrices (%d,%d) (%d,%d)" % ( M,K1,K2,N); return 0
        r = Smatrix(M,N)
        for i in range(M):
            for j in range(N):
                r[ (i,j) ] = strmul( self[(i,0)], other[(0,j)] )
        for k in range(1,K1):
            for i in range(M):
                for j in range(N):
                    r[ (i,j) ]  = stradd( r[(i,j)],
                                          strmul( self[(i,k)], other[(k,j)] ) )
        return r
    def eqns(self,irange,jrange):
        for i in range(min(self.isize,irange)):
            for j in range(min(self.jsize,jrange)):
                yield self[(i,j)]


threesets = [ [ [0,0] ],
              [ [0,0], [0,1] ],
              [ [0,0], [1,0] ],
              [ [0,0], [0,1], [1,0], [1,1] ],
              ]
threediag = [ [ [0,0] ],
              [ [0,0], [0,1], [1,0] ],
              [ [0,0], [0,1], [1,0], [1,1] ],
              ]
threeline = [ [ [0,0] ],
              [ [0,0], [0,1] ],
              ]
threenull = [ [[ 0,0 ]] ]
threecol  = [ [[0,0], [1,0]] ]

head = """
\"\"\"Implements the complete conjugant gradient examples.\"\"\"

from numpy import matrix

from ignition.flame import *
from ignition.utils import flatten_list

# Define the operation
def CG_Op ( %s ):
    return [A * P * D - R * (I - J),
            P * (I - U) - R,
            T(R) * R - W, T(P) * A * P - Z ]

# Define the loop invariant
def CG_Inv ( %s ):
"""
tail1 = """
    known = []

    def unroll_mats(acc, eqn):
        if isinstance(eqn, matrix):
            return acc + flatten_list(eqn.tolist())
        else:
            return acc + [eqn]

    eqns = reduce(unroll_mats, eqns, [])
    return eqns, reduce(unroll_mats, known, [])

# Hold your nose
from ignition.flame.tensors.basic_operators import add_invertible
A = Tensor(\"A\", rank=2)
P_0 = Tensor(\"P_0\", rank=2)
add_invertible(T(P_0) * A * P_0)
add_invertible(T(P_0) * A ** 2 * P_0)


# Define the Partition Objs
"""
tail2 = """
# Generate the algorithm
generate(op=CG_Op, loop_inv=CG_Inv, solver=tensor_solver,
         inv_args=[%s], filetype=\"latex\",
         solution_file=\"alg%d.in\", logic_files=\"cg_logic\")
"""
class Var():
    def __init__(self,name,i,j,prefix,argsrc):
        self.name = name; self.i = i; self.j = j; self.prefix = prefix
        self.argsrc = argsrc
        self.object = Smatrix(i,j)
        if prefix=="Diag_" or (i==1 and j==1):
            self.object.setscalar(self.name)
        elif prefix=="I_":
            self.object.setI(); self.argsrc = "Input"
        elif prefix=="J_":
            self.object.setJ(); self.argsrc = "Input"
        elif prefix=="Upper_":
            self.object.setupper(self.name)
        else:
            self.object.setfull(self.name)
    def line(self):
        if self.argsrc is None:
            argsrc = ""
        else:
            argsrc = ", arg_src=\""+self.argsrc+"\""
        return "%s = iterative_arg(\"%s\",rank=2, part_suffix=\"%s%dx%d\"%s)" % \
               (self.name,self.name,self.prefix,self.i,self.j,argsrc)
    def __repr__(self):
        return self.name+"_"+str(self.i)+"x"+str(self.j)
class PME():
    def __init__(self,**kwargs):
        self.exprs = [] # exprs and invars are coupled: 
        self.invar = [] # expressions and index subsets in them.
        self.vars = []
        self.name = kwargs.pop("name","pme")
        self.workdir = kwargs.pop("workdir",os.getcwd()+"/"+self.name+"-files")
    def addvar(self,name,i,j,prefix="",argsrc="Computed"):
        self.vars.append(Var(name,i,j,prefix,argsrc))
    def varstring(self):
        s = ""
        for v in self.vars: s += v.name+","
        return s
    def append(self,expr,invar=None):
        """Add an Smatrix, probably computed as an expression.
        The optional invar argument describes the subsets that can be
        considered for adding to an invariant."""
        self.exprs.append(expr)
        if invar is None:
            self.invar.append( threesets )
        else: self.invar.append( invar )
    def invariants(self):
        for p in itertools.product( *self.invar ):
            # p is now an array of sets of subscripts:
            # one subset for each consecutive expr
            eqns = []
            for i in range(len(p)):
                for s in self.exprs[i].multisub(p[i]):
                    eqns.append(s)
            yield eqns
    def length(self):
        c = 1
        for p in itertools.product( *self.invar ):
            c += 1
        return c
    def printhead(self,c=0,f=None):
        vs = self.varstring()
        if f is None:
            print head % ( vs, vs )
        else:
            f.write( head % ( vs, vs ) )
    def printtail(self,c=0,f=None):
        vs = self.varstring()
        if f is None:
            print tail1
            for v in self.vars:
                print v.line()
            print tail2 % ( vs, c )
        else:
            f.write( tail1 )
            for v in self.vars:
                f.write( v.line()+"\n" )
            f.write( tail2 % ( vs, c ) )
    def invariants_to_file(self):
        c = 1
        if os.path.exists(self.workdir):
            if os.path.isfile(self.workdir):
                print "Directory name <%s> already taken by file\n" % workdir
                sys.exit(1)
        else:
            os.mkdir(self.workdir)
        cf = open(self.workdir+"/"+self.name+"-commandlines","w")
        for eqns in self.invariants():
            fn = self.invariant_to_file(c,eqns,self.workdir,self.name+"-")
            cf.write("python %s > log-%d 2>&1 &\n" % (fn,c) )
            c += 1
        cf.close()
    def invariant_to_file(self,c,eqns,dir,fh):
        ind = "    "
        print "Invariant",c # ,":",eqns
        fn = fh+str(c)+".py"
        f = open(dir+"/"+fn,'w')
        self.printhead(f=f,c=c)
        for v in self.vars:
            f.write( ind+"%s = %s\n" % ( str(v.object),v.name ) )
        # take apart the objects
        f.write(ind+"eqns = [\n")
        for e in eqns:
            f.write(ind+" "+e+" ,\n")
        f.write(ind+"]\n")
        #
        self.printtail(f=f,c=c)
        f.close()
        return fn
