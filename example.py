

####################################################################################
# In this module, you can compute random (ell,ell)-isogeny from  E0*E0 for fixed p,l.
# p: characteristic
# l: degree
# E0: y^2=x^3+x over F_{p^2}.
# For more details, see "test.py" in the same package.
####################################################################################



if __name__ == "__main__":
    from sage.all import *
    from class_theta_dim1 import Dim1_theta_null
    from func_elliptic    import Ext_mont_pt_to_lv2,Ext_basis_mont_pt_to_lv2
    from func_isogeny     import Codomain_dim2,Evaluate_dim2_special
    from func_E0          import Random_kernel_from_E0E0
    #base field======================================================
    p = 276154505650672190920223
    l=11
    print("\n p:",p)
    print("\n ell:",l)
    Fp2=GF(p**2)
    zeta_4=sqrt(Fp2(-1))
    zeta_8=sqrt(zeta_4)
    sqrt2=sqrt(Fp2(2))
    #elliptic curve===================================================
    E0 = EllipticCurve(Fp2, [0, 0, 0, +1, 0]) #E0:y^2=x^3+x
    print("\n E_0:\n",E0)
    #kernel===================================================
    D= 3**6 * 7 * 11**4 * 19 * 29**3 * 67 * 79
    d=13 * 23 * 37 * 43 * 47**2 * 59 * 61 * 71 * 73 
    [e1_E1,e1_E2],[e2_E1,e2_E2]=Random_kernel_from_E0E0(p,D,d,l,E0)
    print("\n e1=(e1^(1),e1^(2)):\n",(e1_E1,e1_E2))
    print("\n e2=(e2^(1),e2^(2)):\n",(e2_E1,e2_E2))
    e12_E1=e1_E1+e2_E1
    e12_E2=e1_E2+e2_E2
    #point for evaluation===========================================
    #E_1
    ord=19
    x_E1=E0.random_point()
    xpe1_E1=x_E1+e1_E1
    xpe2_E1=x_E1+e2_E1
    #E_2
    x_E2=E0(0)
    xpe1_E2=x_E2+e1_E2
    xpe2_E2=x_E2+e2_E2
    print("\n x=(x^(1),0):\n",(x_E1,x_E2))
    #compute theta coordinate in dim 1 ===================================================
    tc0_E0=[sqrt2,1+zeta_4]
    list_tnp_E1=tc0_E0
    list_tnp_E2=tc0_E0
    #E_1
    tnp_E1=Dim1_theta_null(list_tnp_E1,1,Fp2)
    ext_basis_E1=Ext_basis_mont_pt_to_lv2(list_tnp_E1,e1_E1,e2_E1,p,Fp2)
    ext_x_E1 =Ext_mont_pt_to_lv2(list_tnp_E1,x_E1,e1_E1,e2_E1,p,Fp2)
    #E_1
    tnp_E2=Dim1_theta_null(list_tnp_E2,1,Fp2)
    ext_basis_E2=Ext_basis_mont_pt_to_lv2(list_tnp_E1,e1_E2,e2_E2,p,Fp2)
    ext_x_E2=Ext_mont_pt_to_lv2(list_tnp_E2,x_E2,e1_E2,e2_E2,p,Fp2)
    #isogeny ================================================================================
    #codomain
    tc_f0,(coeff_E1,coeff_E2),(exc_h_lincom_lsq_E1,exc_h_lincom_lsq_E2)=Codomain_dim2(tnp_E1,tnp_E2,ext_basis_E1,ext_basis_E2,l)
    print("\n theta-null point of the codomain:\n", tc_f0.th)
    #evaluation
    tc_fx_s =Evaluate_dim2_special(tnp_E1,ext_basis_E1,ext_x_E1,l,coeff_E1,exc_h_lincom_lsq_E2,coeff_E2)
    print("\n theta coordinate of the image f(x):\n",tc_fx_s.th)
    print("\n")
    
    
    
