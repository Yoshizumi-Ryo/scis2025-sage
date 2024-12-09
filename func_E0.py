

######################################################################
# Functions about E_0: y^2=x^3+x over F_{p^2} where p=3 (mod 4).
# This elliptic curve is supersingular. 
# The j-invariant of E_0 is 1728.
# We will use the endomorpshim ring of E_0.
######################################################################


from sage.all import *
from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius


#Cyclic isogeny======================================================================


def Decomp_degree(N:int):
    """ decompose a given integer to prime factors as list."""
    fac=list(factor(N))
    compo_list=[]
    for i in range(0,len(fac)):
        for counter in range(0,fac[i][1]):
            compo_list.append(fac[i][0])
    return compo_list




def Elliptic_Cyclic(E,ker,P,Q):
    """ 
    Cyclic isogeny from E with kernel "ker".
    Then, compute the image of P,Q.
    """
    deg=order(ker)
    fac_list=Decomp_degree(deg)
    #print("Elliptic isogeny",fac_list)
    s=1
    for i in range(0,len(fac_list)):
        l=fac_list[i]
        k=deg//(s*l)
        decomp_ker=k*ker
        isogeny=EllipticCurveIsogeny(E,decomp_ker)
        E=isogeny.codomain()
        ker=isogeny(ker)
        P=isogeny(P)
        Q=isogeny(Q)
        s*=l
    return E,P,Q


#E_0==============================================================================


def All_4throot_mod(M:int):
    """ compute all integers x such that x^2+1=0 mod M. (0<x<M) """
    assert(M>=2)
    #set of prime factors of M.
    pr_fac_M={list(factor(M))[i][0] for i in range(0,len(list(factor(M))))}
    if M%4==0:
        return set() #empty set.
    for pr in pr_fac_M:
        if pr%4==3:
            return set()  #empty set.
    fac_M=list(factor(M))
    num_pfac=len(fac_M)
    CRT=[fac_M[i][0]**fac_M[i][1] for i in range(0,num_pfac)]
    X_M=[ZZ(Mod(-1,CRT[i]).sqrt())for i in range(0,num_pfac)]
    set_All_4throot_modM=set()
    for bit in subsets(range(0,num_pfac)):
        Y_M=list()
        for i in range(0,num_pfac):
            if i in bit:
                sign=1
            else:
                sign=-1
            Y_M.append((sign*X_M[i])%CRT[i])
        set_All_4throot_modM.add(CRT_list(Y_M,CRT))
    return set_All_4throot_modM

        


def Smoothness(M:int,B:int):
    """ check if M is B-smooth except for one prime factor. """
    prime_under={pr for pr in range(2,B+1) if is_prime(pr)}
    for pr in prime_under:
        while (M%pr==0):
            M=M//pr
    if M==1 or is_prime(M):
        return True
    else:
        return False




def Primitive_Cornacchia(M:int):
    """return r,s such that r^2+s^2=M where gcd(r,s)=1."""
    assert(M>=0)
    if M==0: return True,0,0
    if M==1: return True,1,0
    if M==2: return True,1,1
    if not(Smoothness(M,2**(15))):
        return False,0,0
    set_All_4throot_modM=All_4throot_mod(M)
    for r0 in set_All_4throot_modM:
        if (2*r0)<=M:
            r_0=r0
            assert((r_0**2+1)%M==0)
            r_1=M%r_0
            while (r_1**2>=M):
                r_2=r_0%r_1
                r_0=r_1
                r_1=r_2
            assert(r_1**2 < M)
            if is_square(M-r_1**2):
                s=isqrt(M-r_1**2)
                assert(s**2==(M-r_1**2))
                assert(r_1**2+s**2==M)
                assert(GCD(r_1,s)==1)
                return True,r_1,s
    #Is you come here, you can't find.
    return False,0,0






def FullRepresentInteger(C:int,p:int):
    """ 
    For given integer C, output (x,y,z,t) in B_{p,infty} of norm C.
    Here, (x,y,z,t) means x+y*i+z*(j+k)/2+t*(1+k)/2.
    """
    #assert(C>p)
    assert(p%4==3)
    for counter in range(1,2**(30)):
        cd  = floor(sqrt(4*(C/p)))
        zd  = randint(-cd,cd)
        cdd = floor(sqrt((4*(C/p))-zd**2))
        td  = randint(-cdd,cdd)
        c=4*C-p*(zd**2+td**2)
        TF,xd,yd=Primitive_Cornacchia(c)
        if TF:
            assert(xd**2+yd**2==c)
            if ((xd-td)%2==0) and ((yd-zd)%2==0):
                x=(xd-td)//2
                y=(yd-zd)//2
                z=zd
                t=td
                deg=(xd**2+yd**2+p*(zd**2+td**2))//4
                assert(deg==C)
                return (x,y,z,t)
    #If you come here, you can't find the example.
    assert(False)





def Image_by_RepInt(rep_int:tuple,pt_list:list,end_i):
    """ 
    E_0:y^2=x^3+x / F_p.
    When an endomorphism of E_0 is given as an element of quaternion algebra, 
    compute the image of given points P,Q in E_0.
    Here, rep_int=(x,y,z,t) means x+y*i+z*(j+k)/2+t*(1+k)/2 where j is the Frobenius.
    """
    E0=pt_list[0].curve()
    end_j=EllipticCurveHom_frobenius(E0)
    end_k=end_i*end_j
    [x,y,z,t]=rep_int
    img_list=[]
    for counter in range(0,len(pt_list)):
        P=pt_list[counter]
        hP=P.division_points(2)[0]
        termP_3=(end_i(hP)+end_j(hP))
        termP_4=(hP+end_k(hP))
        img_P=x*P+y*end_i(P)+z*termP_3+t*termP_4
        img_list.append(img_P)
    return img_list




def Pre_Random_Isog_Images(p:int,D:int,d:int,E0):
    """ 
    E0:y^2=x^3+x / F_{p^2}.
    D:smooth integer with D|(p^2+1).
    d: any integer with d<D.
    (P0,Q0) is a basis of E[D]. 
    For d*(D-d)-endomorphism alpha:E_0->E_0,  
    we decompose d-isogeny tau:E_0->E_cd and (D-d)-isogeny E_cd->E_0.
    """
    assert(D>d)
    assert(10*(D-d)*d>p)
    assert((p+1)%D==0)
    assert(GCD(D,d)==1)
    assert(E0.is_supersingular())
    alpha_rep_int=FullRepresentInteger(d*(D-d),p)
    aut_i=E0.isomorphisms(E0)[3]
    P0,Q0=E0.torsion_basis(D)
    [alpha_P0,alpha_Q0]=Image_by_RepInt(alpha_rep_int,[P0,Q0],aut_i)
    assert (alpha_P0 in E0 and alpha_Q0 in E0)
    alpha_P0=alpha_P0*(order(alpha_P0)//D)
    alpha_Q0=alpha_Q0*(order(alpha_Q0)//D)
    return [(D-d)*P0,alpha_P0],[(D-d)*Q0,alpha_Q0]




def Random_kernel_from_E0E0(p:int,D:int,d:int,l:int,E0):
    """ 
    E0:y^2=x^3+x / F_{p^2}.
    D:smooth integer with D|(p^2+1).
    l: prime factor of D.
    Give one basis of random kernel from E0*E0.
    """ 
    assert(D%l==0)
    assert(is_prime(l))
    [kerD_1_1,kerD_1_2],[kerD_2_1,kerD_2_2]=Pre_Random_Isog_Images(p,D,d,E0)
    e1_E1=kerD_1_1*(D//l)
    e1_E2=kerD_1_2*(D//l)
    e2_E1=kerD_2_1*(D//l)
    e2_E2=kerD_2_2*(D//l)
    assert e1_E1.order()==l
    assert e2_E1.order()==l
    assert e1_E2.order()==l
    assert e2_E2.order()==l
    while e1_E1.weil_pairing(e2_E1,l)*e1_E2.weil_pairing(e2_E2,l)!=1:
        e2_E2*=2
    assert e1_E1.weil_pairing(e2_E1,l)*e1_E2.weil_pairing(e2_E2,l)==1
    return [e1_E1,e1_E2],[e2_E1,e2_E2]




