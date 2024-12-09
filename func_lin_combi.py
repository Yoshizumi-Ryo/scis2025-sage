
######################################################################
# In this module, we give some functions for linear combinations. 
# The functions will be used for (l,l)-isogeny calculations.
######################################################################

from sage.all import *
from class_theta_dim1 import Dim1_theta_null


def Set_H_ell(l:int):
    """ define the set H_ell."""
    ld=(l-1)//2
    H_ell={(1,0),(0,1),(1,1)}
    for k2 in range(2,ld+1):
        H_ell.add((0,k2))
    for k2 in range(2,l-1):
        H_ell.add((1,k2))
    for k1 in range(2,ld+1):
        H_ell.add((k1,0))
    for k1 in range(2,l):
        H_ell.add((k1,1))
    for k1 in range(2,l-1):
        H_ell.add((k1,2))
    for k1 in range(2,ld+1):
        for k2 in range(3,l-k1):
            H_ell.add((k1,k2))
    for k1 in range(ld+1,l):
        for k2 in range(3,l-k1+1):
            H_ell.add((k1,k2))
    (len(H_ell)==(l**2-1)//2)
    return H_ell



def Comp_H_ell(l:int):
    """ 
    Complement of H_ell in {(k1,k2) in Z^2 | 0<=k1,k2<l}-{(0,0)}.
    Note that (0,0) is not in comp_H_ell.
    """
    H_ell=Set_H_ell(l)
    c_H_ell=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            if (not (k1,k2) in H_ell) and (k1,k2)!=(0,0):
                if k1==0:
                    assert((0,l-k2) in H_ell)
                    c_H_ell[(0,k2)]=(0,l-k2)
                elif k2==0:
                    assert((l-k1,0) in H_ell)
                    c_H_ell[(k1,0)]=(l-k1,0)
                else:
                    assert((l-k1,l-k2) in H_ell)
                    c_H_ell[(k1,k2)]=(l-k1,l-k2)
    assert (len(c_H_ell)==(l**2-1)//2)
    return c_H_ell




def Half_LinCom_dim1(tc_0:Dim1_theta_null,ext_basis:list,l:int):
    """ 
    Now, ext_basis=[e_1,e_2,e_1+e_2] where e_1,e_2 is a basis of E[l] where E is an elliptic curve.
    Compute the affine theta coordinate of k_1*e_1+k_2*e_2 for (k_1,k_2) in H_l.
    """
    [tc_e1,tc_e2,tc_e12]=ext_basis
    #assert(is_odd(l))
    #assert(is_prime(l))
    ld=(l-1)//2
    #lincom[(k1,k2)]=theta coordinate of (k1*e1+k2*e2).
    lincom={(0,0):tc_0,
            (1,0):tc_e1,
            (0,1):tc_e2,
            (1,1):tc_e12}
    for k2 in range(2,ld+2):
        assert(not (0,k2) in lincom)
        lincom[(0,k2)]=tc_0.Diff_add_dim1(lincom[(0,k2-1)],tc_e2,lincom[(0,k2-2)])
    for k2 in range(2,l-1):
        assert(not (1,k2) in lincom)
        lincom[(1,k2)]=tc_0.Diff_add_dim1(lincom[(1,k2-1)],tc_e2,lincom[(1,k2-2)])
    for k1 in range(2,ld+2):
        assert(not (k1,0) in lincom)
        lincom[(k1,0)]=tc_0.Diff_add_dim1(lincom[(k1-1,0)],tc_e1,lincom[(k1-2,0)])
    for k1 in range(2,l):
        assert(not (k1,1) in lincom)
        lincom[(k1,1)]=tc_0.Diff_add_dim1(lincom[(k1-1,1)],tc_e1,lincom[(k1-2,1)])
    for k1 in range(2,l-1):
        assert(not (k1,2) in lincom)
        lincom[(k1,2)]=tc_0.Diff_add_dim1(lincom[(k1-1,2)],tc_e1,lincom[(k1-2,2)])
    for k1 in range(2,ld+1):
        for k2 in range(3,l-k1):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=tc_0.Diff_add_dim1(lincom[(k1,k2-1)],tc_e2,lincom[(k1,k2-2)])
    for k1 in range(ld+1,l):
        for k2 in range(3,l-k1+1):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=tc_0.Diff_add_dim1(lincom[(k1,k2-1)],tc_e2,lincom[(k1,k2-2)])
    lincom[(ld+1,ld+1)]=tc_0.Diff_add_dim1(lincom[(ld+1,ld)],tc_e2,lincom[(ld+1,ld-1)])
    #assert(len(lincom.keys())==(l**2-1)//2)
    #assert(Set_H_ell(l)==lincom.keys())
    return lincom



    
    
    



def XpLinCom_dim1(tc_0:Dim1_theta_null,ext_basis:list,ext_x:tuple,l:int):
    """ 
    In dim 1, we compute affine theta coordinate of x+k_1*e_1+k_2*e_2 for 0<=k_1,k_2<l.
    Remark that we don't need to use this function if x=0.
    Here, ext_basis=[e_1,e_2,e_1+e_2] and ext_x=[x,x+e_1,x+e_2].
    """
    [tc_e1,tc_e2,tc_e12]=ext_basis
    (tc_x,tc_xpe1,tc_xpe2)=ext_x
    tc_xpe12=tc_0.Three_way_dim1(tc_x,tc_e1,tc_e2,tc_xpe1,tc_e12,tc_xpe2)
    #xplincom[(k1,k2)]=theta coordinate of (x+k1*e1+k2*e2).
    xplincom={}
    xplincom={(0,0):tc_x,
              (1,0):tc_xpe1,
              (0,1):tc_xpe2,
              (1,1):tc_xpe12}
    for k2 in range(2,l+1):
        xplincom[(0,k2)]=tc_0.Diff_add_dim1(xplincom[(0,k2-1)],tc_e2,xplincom[(0,k2-2)])
    for k2 in range(2,l):
        xplincom[(1,k2)]=tc_0.Diff_add_dim1(xplincom[(1,k2-1)],tc_e2,xplincom[(1,k2-2)])
    for k1 in range(2,l+1):
        assert(not (k1,k2) in xplincom.keys())
        xplincom[(k1,0)]=tc_0.Diff_add_dim1(xplincom[(k1-1,0)],tc_e1,xplincom[(k1-2,0)])
    for k2 in range(1,l):
        for k1 in range(2,l):
            assert(not (k1,k2) in xplincom.keys())
            xplincom[(k1,k2)]=tc_0.Diff_add_dim1(xplincom[(k1-1,k2)],tc_e1,xplincom[(k1-2,k2)])
    return xplincom


