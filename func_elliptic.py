

##########################################################################
# functions to transform from Montgomery coordinate to theta coordinate.
##########################################################################

from sage.all import *
from class_theta_dim1 import Dim1_theta,Dim1_theta_null


def Mont_4torsion_to_lv2null(T:list):
    """
    From Montgomery to theta.
    Given 4-torsion point T'_2 s.t. (T'_1;T'_2) is a symeplectic basis with T'_1=(-1:1) in Montgomery coordinate, 
    output theta-null point of the induced theta structure.
    """
    assert(type(T)==list)
    assert(len(T)==2)
    r=T[0] #x-coord.
    s=T[1] #z-coord.
    a=r+s
    b=r-s
    return [a,b]
    


def Mont_pt_to_lv2(lv2tnp:list,P,p:int,field):
    """ 
    From Montgomery coordinate, output level 2 theta coordinate.
    """
    if P==0:
        return Dim1_theta([lv2tnp[0],lv2tnp[1]],1,field)
    a=lv2tnp[0]
    b=lv2tnp[1]
    x=P[0]
    z=P[2]
    tc_0=a*(x-z)
    tc_1=b*(x+z)
    return Dim1_theta([tc_0,tc_1],1,field)


def Ext_basis_mont_pt_to_lv2(list_tnp:list,e1,e2,p:int,field):
    tc_e1 =Mont_pt_to_lv2(list_tnp,e1   ,p,field)
    tc_e2 =Mont_pt_to_lv2(list_tnp,e2   ,p,field)
    tc_e12=Mont_pt_to_lv2(list_tnp,e1+e2,p,field)
    return tc_e1,tc_e2,tc_e12




def Ext_mont_pt_to_lv2(list_tnp:list,x,e1,e2,p:int,field):
    """ 
    For x in A, compute affine lifts of x, x+e_1, x+e_2.
    Here, (e_1,e_2) is a basis of some subgroup ~= (Z/lZ)^2.
    """
    tc_x   =Mont_pt_to_lv2(list_tnp,x   ,p,field)
    tc_xpe1=Mont_pt_to_lv2(list_tnp,x+e1,p,field)
    tc_xpe2=Mont_pt_to_lv2(list_tnp,x+e2,p,field)
    return tc_x,tc_xpe1,tc_xpe2



