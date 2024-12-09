

######################################################################
# In this module, you can compute examples of (ell,ell)-isogeny using.
######################################################################



# #===================================
# import importlib

# import class_theta_dim1
# import class_theta_dim2
# import func_fraction
# import func_elliptic
# import func_isogeny
# import func_E0
# import func_lin_combi

# importlib.reload(class_theta_dim1)
# importlib.reload(class_theta_dim2)
# importlib.reload(func_fraction)
# importlib.reload(func_elliptic)
# importlib.reload(func_isogeny)
# importlib.reload(func_E0)
# importlib.reload(func_lin_combi)
# #===================================



from sage.all import *
from class_theta_dim1 import Dim1_theta_null
from func_elliptic    import Ext_mont_pt_to_lv2,Ext_basis_mont_pt_to_lv2
from func_isogeny     import Codomain_dim2,Evaluate_dim2_general,Evaluate_dim2_general_variant,Evaluate_dim2_special
from func_E0          import Random_kernel_from_E0E0



#base field======================================================


p = 276154505650672190920223
D= 3**6 * 7 * 11**4 * 19 * 29**3 * 67 * 79
d=13 * 23 * 37 * 43 * 47**2 * 59 * 61 * 71 * 73 
l=11

assert(D>d)
assert(is_prime(p))
assert(p%4==3)
assert((p+1)%D==0)
assert(D%l==0)

Fp2.<zeta_4>=GF(p**2,modulus=x**2+1,name="zeta_4")
zeta_8=sqrt(zeta_4)
sqrt2=sqrt(Fp2(2))

#elliptic curve===================================================


#E0:y^2=x^3+x
E0 = EllipticCurve(Fp2, [0, 0, 0, +1, 0])

assert(E0.is_supersingular())
assert(E0.order()==(p+1)**2)
assert(E0.j_invariant()==1728)

tc0_E0=[sqrt2,1+zeta_4]

E1=E0
E2=E0

#kernel====================================


[e1_E1,e1_E2],[e2_E1,e2_E2]=Random_kernel_from_E0E0(p,D,d,l,E0)
            
e12_E1=e1_E1+e2_E1
e12_E2=e1_E2+e2_E2


#point for evaluation===========================================

#E_1
ord=19
x_E1=E1.torsion_basis(ord)[0]
xpe1_E1=x_E1+e1_E1
xpe2_E1=x_E1+e2_E1

#E_2
#x_E2=E2.random_point()
x_E2=E2(0)
xpe1_E2=x_E2+e1_E2
xpe2_E2=x_E2+e2_E2


#compute theta coordinate in dim 1 ===================================================

list_tnp_E1=tc0_E0
list_tnp_E2=tc0_E0

#E_1
tnp_E1=Dim1_theta_null([list_tnp_E1[0],list_tnp_E1[1]],1,Fp2)
ext_basis_E1=Ext_basis_mont_pt_to_lv2(list_tnp_E1,e1_E1,e2_E1,p,Fp2)
ext_x_E1 =Ext_mont_pt_to_lv2(list_tnp_E1,x_E1,e1_E1,e2_E1,p,Fp2)


#E_2
tnp_E2=Dim1_theta_null([list_tnp_E2[0],list_tnp_E2[1]],1,Fp2)
ext_basis_E2=Ext_basis_mont_pt_to_lv2(list_tnp_E2,e1_E2,e2_E2,p,Fp2)
ext_x_E2=Ext_mont_pt_to_lv2(list_tnp_E2,x_E2,e1_E2,e2_E2,p,Fp2)

#codomain of isogeny=========================================================


tc_f0,(coeff_E1,coeff_E2),(exc_h_lincom_lsq_E1,exc_h_lincom_lsq_E2)=Codomain_dim2(tnp_E1,tnp_E2,ext_basis_E1,ext_basis_E2,l)


#evaluation (general)========================================================================

tc_fx_g=Evaluate_dim2_general(tnp_E1,ext_basis_E1,ext_x_E1,coeff_E1,
                                tnp_E2,ext_basis_E2,ext_x_E2,coeff_E2,
                                l)
                              
assert tc_f0.Is_order(tc_fx_g,ord)

#evaluation (general variant )========================================================================

tc_fx_g_v=Evaluate_dim2_general_variant(tnp_E1,ext_basis_E1,ext_x_E1,coeff_E1,
                                tnp_E2,ext_basis_E2,ext_x_E2,coeff_E2,
                                l)
                              
assert tc_f0.Is_order(tc_fx_g_v,ord)

#evaluation (special)========================================================================                   
                    
                    
tc_fx_s =Evaluate_dim2_special(tnp_E1,ext_basis_E1,ext_x_E1,l,coeff_E1,exc_h_lincom_lsq_E2,coeff_E2)

assert tc_f0.Is_order(tc_fx_s,ord)

#check ================================================================================================

assert tc_fx_g.Peq(tc_fx_g_v)
assert tc_fx_g_v.Peq(tc_fx_s)
assert tc_fx_g.Peq(tc_fx_s)

#====================================================================================================
#====================================================================================================