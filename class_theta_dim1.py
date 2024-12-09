
######################################################################
# In this module, we define some class about theta coordinate.
######################################################################

from sage.all import *
from class_theta_dim2 import Coord



class Dim1_theta:
    """ 
    theta coordinate in dimension 1.
    we hold the coordinate by 3 elements [n0,n1,d] s.t. the i-th component is n_i/d.
    """
    def __init__(self,numer:list,denom,base_field):
        assert(len(numer)==2)
        assert(denom!=0) #denominator is not 0.
        self.ch    =base_field.characteristic()
        self.field =base_field
        self.th0   =base_field(numer[0])
        self.th1   =base_field(numer[1])
        self.numer =[self.th0,self.th1]
        self.denom =base_field(denom)
        self.th    =[self.th0,self.th1,self.denom]
        
        
    def Num_square(self):
        """ compute squares of only numerators. this require 2-times squares on the base field. """ 
        try:
            return self.num_sq
        except AttributeError:
            self.num_sq=[self.th0**2,self.th1**2]
            return self.num_sq


    def Product_theta(self,tc):
        """ 
        compute 2-dimensional theta coordinate from 2 1-dimensional theta. 
        this require 5-times multiplications on the base field
        """
        numer=[self.th0*tc.th0, self.th1*tc.th0, self.th0*tc.th1, self.th1*tc.th1]
        denom=self.denom*tc.denom
        return Coord(numer,denom,self.field)
        
  
    def Aeq(self,tc):
        """ equal as affine theta coordinate."""
        if self.th0/self.denom==tc.th0/tc.denom and self.th1/self.denom==tc.th1/tc.denom:
            return True 
        else:
            return False  

    def Peq(self,tc):
        """ equal as projective theta coordinate."""
        if self.th0/self.th1==tc.th0/tc.th1:
            return True 
        else:
            return False  



######################################################################
######################################################################



class Dim1_theta_null(Dim1_theta):
    """ 
    class for theta-null point in dimension 1.
    We hold the coordinate by 2 elements a,b. However, we input as [a,b,1].
    """
    @staticmethod
    def Common_denominator(frac_1,frac_2):
        """ compute common denominator. i,e, given a/b and c/d, output (a*d)/(b*d) and (d*c)/(b*d). """
        assert(len(frac_1)==2)
        assert(len(frac_2)==2)
        common_denom=frac_1[1]*frac_2[1]
        return [[frac_1[0]*frac_2[1],common_denom],[frac_1[1]*frac_2[0],common_denom]]


    def Dual_sq(self):
        """ for the theta-null point (a,b), output a^2+b^2, a^2-b^2."""
        try:
            return self.dual_sq
        except AttributeError:
            a_sq=self.th0**2
            b_sq=self.th1**2
            A_sq=a_sq+b_sq
            B_sq=a_sq-b_sq
            self.dual_sq=[A_sq,B_sq,A_sq*B_sq]
            return self.dual_sq
        


    def Diff_add_dim1(self,tc_x:Dim1_theta,tc_y:Dim1_theta,tc_xmy:Dim1_theta):
        """ differential addition in dimension 1."""
        [Asq,Bsq,AsqBsq]=self.Dual_sq()
        tc_x_sq=tc_x.Num_square()
        tc_y_sq=tc_y.Num_square()
        dxdy_sq=(tc_x.denom*tc_y.denom)**2
        F_0=[(tc_x_sq[0]+tc_x_sq[1])*(tc_y_sq[0]+tc_y_sq[1]),Asq]
        F_1=[(tc_x_sq[0]-tc_x_sq[1])*(tc_y_sq[0]-tc_y_sq[1]),Bsq]
        n_0=F_0[0]*F_1[1]
        n_1=F_0[1]*F_1[0]
        d=AsqBsq
        G_0=[(n_0+n_1)*tc_xmy.denom,tc_xmy.th0]
        G_1=[(n_0-n_1)*tc_xmy.denom,tc_xmy.th1]
        [[alpha_xpy,d_xpy],[beta_xpy,d_xpy]]=Dim1_theta_null.Common_denominator(G_0,G_1)
        d_xpy*=(2*d*dxdy_sq)
        return Dim1_theta([alpha_xpy,beta_xpy],d_xpy,self.field)
    


    
    def Three_way_dim1(self,tc_x:Dim1_theta,tc_y:Dim1_theta,tc_z:Dim1_theta,
                       tc_xpy:Dim1_theta,tc_ypz:Dim1_theta,tc_zpx:Dim1_theta):
        """ three-way addition in dimension 1."""
        a=self.th0
        b=self.th1
        alpha_y_alpha_z     =tc_y.th0*tc_z.th0
        a_alpha_ypz         =a*tc_ypz.th0
        alpha_zpx_alpha_xpy =tc_zpx.th0*tc_xpy.th0
        beta_y_beta_z       =tc_y.th1*tc_z.th1
        b_beta_ypz          =b*tc_ypz.th1
        beta_zpx_beta_xpy   =tc_zpx.th1*tc_xpy.th1 
        ddd_1=tc_xpy.denom*tc_ypz.denom*tc_zpx.denom
        ddd_2=tc_x.denom*tc_y.denom*tc_z.denom
        F_0=[(a_alpha_ypz+b_beta_ypz)*(alpha_zpx_alpha_xpy+beta_zpx_beta_xpy),(alpha_y_alpha_z+beta_y_beta_z)]
        F_1=[(a_alpha_ypz-b_beta_ypz)*(alpha_zpx_alpha_xpy-beta_zpx_beta_xpy),(alpha_y_alpha_z-beta_y_beta_z)]
        [[n_0,d],[n_1,d]]=Dim1_theta_null.Common_denominator(F_0,F_1)
        G_0=[(n_0+n_1)*ddd_2,2*d*tc_x.th0]
        G_1=[(n_0-n_1)*ddd_2,2*d*tc_x.th1]
        [[alpha_xyz,d_xyz],[beta_xyz,d_xyz]]=Dim1_theta_null.Common_denominator(G_0,G_1)
        d_xyz*=ddd_1
        return Dim1_theta([alpha_xyz,beta_xyz],d_xyz,self.field)
    
  
    def Mult_dim1(self,tc_x:Dim1_theta,n:int):
        """ Not optimized!!!! """
        tc_mx=dict() #tc_mx[m] is theta of m*x.
        tc_mx[0]=self 
        tc_mx[1]=tc_x 
        for m in range(2,n+1):
            tc_mx[m]=self.Diff_add_dim1(tc_mx[m-1],tc_x,tc_mx[m-2])
        return tc_mx[n]


    def Is_order_dim1(self,tc_x:Dim1_theta,n:int):
        """ 
        Not optimized!!!! 
        For a given integer n, we check if x is of order n.
        """
        assert n>=1
        return self.Peq(self.Mult_dim1(tc_x,n))




######################################################################
######################################################################
