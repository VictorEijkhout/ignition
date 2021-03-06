#!/usr/bin/env python

# straight CG

from pme import *

if __name__== "__main__":
    pme = PME(name="cg-ortho")
    A = Smatrix(1,1); A.setscalar("A")
    pme.addvar("A",1,1,argsrc="Input")
    R = Smatrix(1,3); R.setfull("R")
    pme.addvar("R",1,3,argsrc="Overwrite")
    P = Smatrix(1,3); P.setfull("P")
    pme.addvar("P",1,3)
    D = Smatrix(3,3); D.setscalar("D")
    pme.addvar("D",3,3,prefix="Diag_")
    I = Smatrix(3,3); I.setI()
    pme.addvar("I",3,3,prefix="I_")
    J = Smatrix(3,3); J.setJ()
    pme.addvar("J",3,3,prefix="J_")
    pme.append( R*J-R*I-A*P*D, threenull )

    U = Smatrix(3,3); U.setunitupper("U")
    pme.addvar("U",3,3,"Upper_")
    pme.append( P*U-R, threeline )

    W = Smatrix(3,3); W.setscalar("W")
    pme.addvar("W",3,3,"Diag_")
    pme.append( R.t()*R-W, threediag )

    Z = Smatrix(3,3); Z.setscalar("Z")
    pme.addvar("Z",3,3,"Diag_")
    pme.append( P.t()*A*P-Z, threediag )

    print "\nGoing to try",pme.length(),"sets"
    pme.invariants_to_file()
