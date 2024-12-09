
######################################################################################################################
# In this module, we give some functions for (ell,ell)-isogeny calculation where ell is odd prime number.
# There are 3 types functions: Codomain, Evaluation(general), Evaluation(special). 
######################################################################################################################

from sage.all import *
from func_lin_combi   import Set_H_ell,Half_LinCom_dim1,XpLinCom_dim1,Comp_H_ell
from func_fraction    import Dict_common_denom_len5
from class_theta_dim1 import Dim1_theta,Dim1_theta_null
from class_theta_dim2 import Coord,NullCoord

           
def Prod_nnd_nnd(nnd1:list,nnd2:list):
    """
    "nnd" is 1-dimensional theta coordinate holded as [numerator at 0,numerator at 1,(common) denominator].
    For nnd1, nnd2, we give the product theta coordinate in dim 2.
    """
    assert len(nnd1)==3
    assert len(nnd2)==3
    prod=[nnd1[0]*nnd2[0],
          nnd1[1]*nnd2[0],
          nnd1[0]*nnd2[1],
          nnd1[1]*nnd2[1],
          nnd1[2]*nnd2[2]]
    return prod
    
##########################################################################################

def Lmd_lsq_dim1(tc_ld_e:Dim1_theta,tc_ldp1_e:Dim1_theta):
    """ compute lmd^{\l} to use normalization of affine lift of basis."""
    assert (tc_ld_e.Peq(tc_ldp1_e))
    lmd_lsq_1=[tc_ld_e.th0*tc_ldp1_e.denom,tc_ld_e.denom*tc_ldp1_e.th0]
    lmd_lsq_2=[tc_ld_e.th1*tc_ldp1_e.denom,tc_ld_e.denom*tc_ldp1_e.th1]
    assert lmd_lsq_1[0]/lmd_lsq_1[1]==lmd_lsq_2[0]/lmd_lsq_2[1]
    return lmd_lsq_1




def Coefficient_codomain_dim1(ld_e1 :Dim1_theta,ldp1_e1 :Dim1_theta,
                                      ld_e2 :Dim1_theta,ldp1_e2 :Dim1_theta,
                                      ld_e12:Dim1_theta,ldp1_e12:Dim1_theta,
                                      l:int):
    """ compute coefficients to use the normalization for codomain calculation in dimension 1."""
    #pick up.
    lmd1_lsq =Lmd_lsq_dim1(ld_e1 ,ldp1_e1)
    lmd2_lsq =Lmd_lsq_dim1(ld_e2 ,ldp1_e2)
    lmd12_lsq=Lmd_lsq_dim1(ld_e12,ldp1_e12)
    lmd_div_lsq=[lmd12_lsq[0]*lmd1_lsq[1]*lmd2_lsq[1],lmd12_lsq[1]*lmd1_lsq[0]*lmd2_lsq[0]]
    #power.
    lmd1_lsq_pow   =[[lmd1_lsq[u]**(m**2) for u in range(0,2)] for m in range(0,l)]
    lmd2_lsq_pow   =[[lmd2_lsq[u]**(m**2) for u in range(0,2)] for m in range(0,l)]
    lmd_div_lsq_pow=[[[lmd_div_lsq[u]**(m1*m2) for u in range(0,2)] for m2 in range(0,l)]for m1 in range(0,l)]
    #product.
    coeff=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff[(k1,k2)]=[lmd1_lsq_pow[k1][u]*lmd2_lsq_pow[k2][u]*lmd_div_lsq_pow[k1][k2][u] for u in range(0,2)]
    return coeff






def Normalization_codomain_dim1(h_lincom:dict,coeff:dict,l:int):
    """ normalization in codomain in dimension 1."""
    #take l-th power.
    h_lincom_lsq=dict()
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        h_lincom_lsq[(k1,k2)]=[h_lincom[(k1,k2)].th0**l, 
                               h_lincom[(k1,k2)].th1**l, 
                               h_lincom[(k1,k2)].denom**l]
    #normalization.
    excellent_h_lincom_lsq=dict()
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        excellent_h_lincom_lsq[(k1,k2)]=[h_lincom_lsq[(k1,k2)][0]*coeff[(k1,k2)][0],
                                         h_lincom_lsq[(k1,k2)][1]*coeff[(k1,k2)][0],
                                         h_lincom_lsq[(k1,k2)][2]*coeff[(k1,k2)][1]] 
    return excellent_h_lincom_lsq




def Codomain_dim2(tc_0_E1:Dim1_theta_null,tc_0_E2:Dim1_theta_null,ext_basis_E1:list,ext_basis_E2:list,l:int):
    """ 
    Compute theta-null point of the codomain for (l,l)-isongey.
    Here, e1,e2 is a basis of the kernel in (E1*E2)[l].
    tc_0_E1 is the theta-null point of E1 and ext_basis_E1 is the first component of e1,e2,e1+e2.
    """
    h_lincom_E1=Half_LinCom_dim1(tc_0_E1,ext_basis_E1,l)
    h_lincom_E2=Half_LinCom_dim1(tc_0_E2,ext_basis_E2,l)
    ld=(l-1)//2
    coeff_E1=Coefficient_codomain_dim1(h_lincom_E1[(ld,0)] ,h_lincom_E1[(ld+1,0)],
                                               h_lincom_E1[(0,ld)] ,h_lincom_E1[(0,ld+1)],
                                               h_lincom_E1[(ld,ld)],h_lincom_E1[(ld+1,ld+1)],
                                               l)
    coeff_E2=Coefficient_codomain_dim1(h_lincom_E2[(ld,0)] ,h_lincom_E2[(ld+1,0)],
                                               h_lincom_E2[(0,ld)] ,h_lincom_E2[(0,ld+1)],
                                               h_lincom_E2[(ld,ld)],h_lincom_E2[(ld+1,ld+1)],
                                               l)
    exc_h_lincom_lsq_E1=Normalization_codomain_dim1(h_lincom_E1,coeff_E1,l)
    exc_h_lincom_lsq_E2=Normalization_codomain_dim1(h_lincom_E2,coeff_E2,l)
    exc_h_lincom_lsq_E12=dict()
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        ex_E1=exc_h_lincom_lsq_E1[(k1,k2)]
        ex_E2=exc_h_lincom_lsq_E2[(k1,k2)]
        exc_h_lincom_lsq_E12[(k1,k2)]=Prod_nnd_nnd(ex_E1,ex_E2)
    exc_h_lincom_lsq_E12=Dict_common_denom_len5(exc_h_lincom_lsq_E12)[0]
    tc_f0=[exc_h_lincom_lsq_E12[(0,0)][i] for i in range(0,4)]
    for (k1,k2) in Set_H_ell(l):
        for i in range(0,4):
            tc_f0[i]+=2*exc_h_lincom_lsq_E12[(k1,k2)][i]
    return NullCoord(tc_f0,1,tc_0_E1.field),(coeff_E1,coeff_E2),(exc_h_lincom_lsq_E1,exc_h_lincom_lsq_E2)





####################################################################################################################################
####################################################################################################################################


def Quasi_Mudivlmd_lsq_dim1(tc_x:Dim1_theta,tc_xple:Dim1_theta):
    """ compute (mu/lmd)^l*(lmd^l)^l. """
    value_0=[tc_x.th0*tc_xple.denom,tc_x.denom*tc_xple.th0]
    value_1=[tc_x.th1*tc_xple.denom,tc_x.denom*tc_xple.th1]
    assert value_0[0]/value_0[1]==value_1[0]/value_1[1]
    return value_0



def Coefficient_evaluate_dim2(coeff_E1:dict,tc_x_E1:Dim1_theta,tc_xple1_E1:Dim1_theta,tc_xple2_E1:Dim1_theta,
                                              coeff_E2:dict,tc_x_E2:Dim1_theta,tc_xple1_E2:Dim1_theta,tc_xple2_E2:Dim1_theta,
                                              l:int):
    """ 
    coefficient appearing in the formula to compute evaluation. 
    we reuse "coeff" we computed in codomain.
    """
    lmd1_lsq=[coeff_E1[(1,0)][u]*coeff_E2[(1,0)][u] for u in range(0,2)]
    lmd2_lsq=[coeff_E1[(0,1)][u]*coeff_E2[(0,1)][u] for u in range(0,2)]
    #E1
    qmdl1_E1=Quasi_Mudivlmd_lsq_dim1(tc_x_E1,tc_xple1_E1)
    qmdl2_E1=Quasi_Mudivlmd_lsq_dim1(tc_x_E1,tc_xple2_E1)
    #E2
    qmdl1_E2=Quasi_Mudivlmd_lsq_dim1(tc_x_E2,tc_xple1_E2)
    qmdl2_E2=Quasi_Mudivlmd_lsq_dim1(tc_x_E2,tc_xple2_E2)
    #E1*E2
    mu1_div_lmd1=[qmdl1_E1[0]*qmdl1_E2[0]*lmd1_lsq[1]**l,qmdl1_E1[1]*qmdl1_E2[1]*lmd1_lsq[0]**l]
    mu2_div_lmd2=[qmdl2_E1[0]*qmdl2_E2[0]*lmd2_lsq[1]**l,qmdl2_E1[1]*qmdl2_E2[1]*lmd2_lsq[0]**l]
    #E1*E2
    mu1_div_lmd1=[qmdl1_E1[u]*qmdl1_E2[u]*lmd1_lsq[(u+1)%2]**l for u in range(0,2)]
    mu2_div_lmd2=[qmdl2_E1[u]*qmdl2_E2[u]*lmd2_lsq[(u+1)%2]**l for u in range(0,2)]
    #powers
    mu1_div_lmd1_pow=[[mu1_div_lmd1[u]**m1 for u in range(0,2)] for m1 in range(0,l)]
    mu2_div_lmd2_pow=[[mu2_div_lmd2[u]**m2 for u in range(0,2)] for m2 in range(0,l)]
    eva_coeff=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            eva_coeff[(k1,k2)]=[coeff_E1[(k1,k2)][u]*coeff_E2[(k1,k2)][u]*mu1_div_lmd1_pow[k1][u]*mu2_div_lmd2_pow[k2][u] for u in range(0,2)]
    return eva_coeff






def Normalization_evaluation_dim2(xplincom_E1:dict,xplincom_E2:dict,eva_coeff:dict,l:int):
    """ 
    normalization in evaluation in dimension 2.
    """
    xplincom_E12_lsq=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            xplincom_E12_k1k2=xplincom_E1[(k1,k2)].Product_theta(xplincom_E2[(k1,k2)])
            xplincom_E12_lsq[(k1,k2)]=[xplincom_E12_k1k2.th[i]**l for i in range(0,5)]
    excellent_xplincom_lsq=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            excellent_xplincom_lsq[(k1,k2)]=[xplincom_E12_lsq[(k1,k2)][0]*eva_coeff[(k1,k2)][0],
                                             xplincom_E12_lsq[(k1,k2)][1]*eva_coeff[(k1,k2)][0],
                                             xplincom_E12_lsq[(k1,k2)][2]*eva_coeff[(k1,k2)][0],
                                             xplincom_E12_lsq[(k1,k2)][3]*eva_coeff[(k1,k2)][0],
                                             xplincom_E12_lsq[(k1,k2)][4]*eva_coeff[(k1,k2)][1]] 
    return excellent_xplincom_lsq




def Evaluate_dim2_general(tc_0_E1:Dim1_theta_null,ext_basis_E1:list,ext_x_E1:tuple,coeff_E1:dict,
                          tc_0_E2:Dim1_theta_null,ext_basis_E2:list,ext_x_E2:tuple,coeff_E2:dict,
                          l:int):
    """ 
    evaluation for the general case, i.e., x=(x_E1,x_E2).
    The input "coeff" is reused value we computed when we calculated codomain.
    """
    #linear combination in dim 1 ---------------------------------------------------------------
    xplincom_E1=XpLinCom_dim1(tc_0_E1,ext_basis_E1,ext_x_E1,l)
    xplincom_E2=XpLinCom_dim1(tc_0_E2,ext_basis_E2,ext_x_E2,l)
    #calculate coefficient in dim 2 ---------------------------------------------------------------
    eva_coeff=Coefficient_evaluate_dim2(coeff_E1,ext_x_E1[0],xplincom_E1[(l,0)],xplincom_E1[(0,l)],
                                                        coeff_E2,ext_x_E2[0],xplincom_E2[(l,0)],xplincom_E2[(0,l)],
                                                        l)
    #compute excellent lifts in dim2 ---------------------------------------------------------------
    excellent_xplincom_lsq=Normalization_evaluation_dim2(xplincom_E1,xplincom_E2,eva_coeff,l)
    #take common denominator ---------------------------------------------------------------
    excellent_xplincom_lsq=Dict_common_denom_len5(excellent_xplincom_lsq)[0]
    #compute f(x)---------------------------------------------------------------
    tc_fx=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            for i in range(0,4):
                tc_fx[i]+=excellent_xplincom_lsq[(k1,k2)][i] 
    return Coord(tc_fx,1,tc_0_E1.field)




def Evaluate_dim2_general_variant(tc_0_E1:Dim1_theta_null,ext_basis_E1:list,ext_x_E1:tuple,coeff_E1:dict,
                            tc_0_E2:Dim1_theta_null,ext_basis_E2:list,ext_x_E2:tuple,coeff_E2:dict,
                            l:int):
    #linear combination in dim 1 ---------------------------------------------------------------
    xplincom_E1=XpLinCom_dim1(tc_0_E1,ext_basis_E1,ext_x_E1,l)
    xplincom_E2=XpLinCom_dim1(tc_0_E2,ext_basis_E2,ext_x_E2,l)
    #calculate coefficient in dim 1 ---------------------------------------------------------------
    eva_coeff_E1=Coefficient_evaluate_dim1(l,coeff_E1,ext_x_E1[0],xplincom_E1[(l,0)],xplincom_E1[(0,l)])
    eva_coeff_E2=Coefficient_evaluate_dim1(l,coeff_E2,ext_x_E2[0],xplincom_E2[(l,0)],xplincom_E2[(0,l)])
    #compute excellent lifts in dim1---------------------------------------------------------------
    exc_xplincom_lsq_E1=Normalization_evaluation_dim1(xplincom_E1,eva_coeff_E1,l)
    exc_xplincom_lsq_E2=Normalization_evaluation_dim1(xplincom_E2,eva_coeff_E2,l)
    #take product on E1*E2 ---------------------------------------------------------------
    exc_xplincom_lsq_E12=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            exc_xplincom_lsq_E12[(k1,k2)]=Prod_nnd_nnd(exc_xplincom_lsq_E1[(k1,k2)],exc_xplincom_lsq_E2[(k1,k2)])
    assert(len(exc_xplincom_lsq_E12)==l**2)
    #take common denominator ---------------------------------------------------------------
    exc_xplincom_lsq_E12=Dict_common_denom_len5(exc_xplincom_lsq_E12)[0]
    #compute f(x)---------------------------------------------------------------
    tc_fx=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            for i in range(0,4):
                tc_fx[i]+=exc_xplincom_lsq_E12[(k1,k2)][i] 
    return Coord(tc_fx,1,tc_0_E1.field)






####################################################################################################################################
####################################################################################################################################




def Coefficient_evaluate_dim1(l:int,coeff:dict,tc_x:Dim1_theta,tc_xple1:Dim1_theta,tc_xple2:Dim1_theta):
    """ 
    this function is used in the special codomain to compute coefficient of the non-zero component of x.
    For example, if x=(x_E1,0_E2), this function is used to normalize x_E1.
    """
    lmd1_lsq=coeff[(1,0)]
    lmd2_lsq=coeff[(0,1)]    
    qmdl1=Quasi_Mudivlmd_lsq_dim1(tc_x,tc_xple1)
    qmdl2=Quasi_Mudivlmd_lsq_dim1(tc_x,tc_xple2)
    #modify
    mdl1=[qmdl1[u]*(lmd1_lsq[(u+1)%2]**l) for u in range(0,2)]
    mdl2=[qmdl2[u]*(lmd2_lsq[(u+1)%2]**l) for u in range(0,2)]
    #powers
    mdl1_pow=[[mdl1[u]**m for u in range(0,2)] for m in range(0,l)]
    mdl2_pow=[[mdl2[u]**m for u in range(0,2)] for m in range(0,l)]
    eva_coeff=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            eva_coeff[(k1,k2)]=[coeff[(k1,k2)][u]*mdl1_pow[k1][u]*mdl2_pow[k2][u] for u in range(0,2)]
    return eva_coeff
  




  


def Normalization_evaluation_dim1(xplincom:dict,eva_coeff:dict,l:int):
    """ 
    normalization in evaluation in dimension 1.
    this function is used in the special evaluation for the non-zero component.
    """
    exc_xplincom_lsq=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            exc_xplincom_lsq[(k1,k2)]=[(xplincom[(k1,k2)].th0**l)  *eva_coeff[(k1,k2)][0],
                                       (xplincom[(k1,k2)].th1**l)  *eva_coeff[(k1,k2)][0],
                                       (xplincom[(k1,k2)].denom**l)*eva_coeff[(k1,k2)][1]]
    return exc_xplincom_lsq


  




def Evaluate_dim2_special(tc_0_E1:Dim1_theta_null,ext_basis_E1:list,ext_x_E1:tuple,l:int,coeff_E1:dict,exc_h_lincom_lsq_E2:dict,coeff_E2:dict):
    """ 
    Evaluation for the case of x=(x_E1,0_E2).
    """
    #linear combination on E1. ---------------------------------------------------------------
    xplincom_E1=XpLinCom_dim1(tc_0_E1,ext_basis_E1,ext_x_E1,l)
    #calculate coefficient on E1. -------------------------------------------------------------------
    eva_coeff_E1=Coefficient_evaluate_dim1(l,coeff_E1,ext_x_E1[0],xplincom_E1[(l,0)],xplincom_E1[(0,l)])
    #compute excellent lifts on E1. -------------------------------------------------------------
    exc_xplincom_lsq_E1=Normalization_evaluation_dim1(xplincom_E1,eva_coeff_E1,l)
    #take product on E1*E2. -------------------------------------------------------------------------------
    exc_xplincom_lsq_E12=dict()
    for (k1,k2) in Set_H_ell(l)|{(0,0)}:
        exc_xplincom_lsq_E12[(k1,k2)]=Prod_nnd_nnd(exc_xplincom_lsq_E1[(k1,k2)],exc_h_lincom_lsq_E2[(k1,k2)])
    c_H_ell=Comp_H_ell(l)
    for (k1,k2) in c_H_ell: 
        (k1d,k2d)=c_H_ell[(k1,k2)]
        exc_xplincom_lsq_E12[(k1,k2)]=Prod_nnd_nnd(exc_xplincom_lsq_E1[(k1,k2)],exc_h_lincom_lsq_E2[(k1d,k2d)])
    assert(len(exc_xplincom_lsq_E12)==l**2)
    #take common denominator --------------------------------------------------------------------------------
    exc_xplincom_lsq_E12=Dict_common_denom_len5(exc_xplincom_lsq_E12)[0]
    #compute f(x)------------------------------------------------------------------------------------
    tc_fx=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            for i in range(0,4):
                tc_fx[i]+=exc_xplincom_lsq_E12[(k1,k2)][i] 
    return Coord(tc_fx,1,tc_0_E1.field)

