

##################################################################################
# Useful functions to treat fractions.
##################################################################################

from sage.all import *


def Binary_Exp(N:int,i:int):
    """ 
    For given integers N,i, 
    give binary representation of N with i-digits.
    """
    N_list=ZZ(N).digits(2)
    N_tup=tuple(N_list)
    assert(len(N_tup)<=i)
    for counter in range(0,i-len(N_tup)):
        N_tup=N_tup+tuple([0])
    return N_tup
        
    

def Common_Denom(denoms:list):
    assert(all(denoms[i]!=0 for i in range(0,len(denoms))))
    N=len(denoms)    
    d_list=ZZ(N-1).digits(2)
    n=len(d_list)
    a={}
    b={}
    a={Binary_Exp(m,n):denoms[m] for m in range(0,N)}
    def B_val(t:tuple):
        return sum([t[j]*(2**j) for j in range(0,len(t))])
    for i in range(1,n+1):
        d_i_list=tuple([d_list[j] for j in range(i,n)])
        assert(len(d_i_list)==n-i)
        d_i_val=B_val(d_i_list)
        #(d_i_list==Binary_Exp(d_i_val,n-i))
        for e_val in range(0,d_i_val+1):
            e_list=Binary_Exp(e_val,n-i)
            if e_list==d_i_list and d_list[i-1]==0:
                a[e_list]=a[tuple([0])+e_list]
            else:
                a[e_list]=a[tuple([0])+e_list]*a[tuple([1])+e_list]
    total_product=a[()]
    b[tuple([0])]=a[tuple([1])]
    b[tuple([1])]=a[tuple([0])]
    for i in reversed(range(0,n-1)):
        d_i_list=tuple([d_list[j] for j in range(i,n)])
        d_i_val=B_val(d_i_list)
        for e_val in range(0,d_i_val+1):
            e_list=Binary_Exp(e_val,n-i)
            if e_list==d_i_list and e_list[0]==0 and d_i_list[0]==0:
                r_e_list=tuple([e_list[i] for i in range(1,len(e_list))])
                b[e_list]=b[r_e_list]
            else:
                r_e_list=tuple([e_list[i] for i in range(1,len(e_list))])
                c=(e_list[0]+1)%2
                c_e_list=tuple([c])+r_e_list
                b[e_list]=b[r_e_list]*a[c_e_list]
    b_list=[b[Binary_Exp(i,n)] for i in range(0,N)]
    assert(total_product!=0)
    return b_list,total_product




def Common_denom_frac(fracs:list):
    assert(all(len(fracs[i])==2 for i in range(0,len(fracs))))
    assert(all(fracs[i][1]!=0   for i in range(0,len(fracs))))
    denoms=[fracs[i][1] for i in range(0,len(fracs))]
    numes,den=Common_Denom(denoms)
    new_fracs=[[(numes[i]*fracs[i][0]),den] for i in range(0,len(fracs))]
    assert(den!=0)
    return new_fracs,den




# def Dict_common_denom(lincom:dict):
#     assert (all(type(lincom[key])==Coord for key in lincom))
#     list_lincom=sorted(lincom.items())
#     k=len(list_lincom)
#     list_val=[(list_lincom[i][1]).denom for i in range(0,k)]
#     list_num,den=Common_Denom(list_val)
#     dict_num={list_lincom[i][0]:list_num[i] for i in range(0,k)}
#     for k in lincom.keys():
#         lincom[k].numer=[dict_num[k]*lincom[k].numer[j] for j in range(0,4)]
#         lincom[k].denom=den
#     assert(den!=0)
#     return lincom,den



def Dict_common_denom_len5(lincom:dict):
    """ 
    the type of input "lincom" is "dict".
    for all keys of the dict, lincom[key]=[n0,n1,n2,n3,d] of length 5. 
    here, "n" means "numerator", "d" means "denominator".
    We output new_lincom:dict s.t. new_lincom[key][i]/new_lincom[key][4]==lincom[key][i]/lincom[key][4] for i=0,1,2,3.
    in addition, for all key, new_lincom[key][4] is the same.
    """
    assert (all(len(lincom[key])==5 for key in lincom))
    list_lincom=sorted(lincom.items())
    k=len(list_lincom)
    list_denoms=[list_lincom[i][1][4] for i in range(0,k)]
    list_mult_num,den=Common_Denom(list_denoms)
    dict_mult_num={list_lincom[i][0]:list_mult_num[i] for i in range(0,k)}
    for key in lincom.keys():
        lincom[key]=[dict_mult_num[key]*lincom[key][i] for i in range(0,4)]
        lincom[key].append(den)
    return lincom,den




def Frac_add(frac1:list,frac2:list):
    """ 
    Addition of fractions.
    For given two fractions a/b, c/d, 
    give a/b+c/d=(ad+bc)/bd.
    """
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    return [frac1[0]*frac2[1]+frac1[1]*frac2[0],frac1[1]*frac2[1]]




def Frac_sub(frac1:list,frac2:list):
    """ 
    Subtraction of fractions.
    For given two fractions a/b, c/d, 
    give a/b-c/d=(ad-bc)/bd.
    """
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    if frac1[0]==0:
        return [-frac2[0],frac2[1]]
    elif frac2[0]==0:
        return frac1
    elif frac1[1]==frac2[1]:
        return [frac1[0]-frac2[0],frac1[1]]
    else:
        return [frac1[0]*frac2[1]-frac1[1]*frac2[0],frac1[1]*frac2[1]]



def Frac_mult(frac1:list,frac2:list):
    """ 
    Multiplication of fractions.
    For given two fractions a/b, c/d, 
    give a/b*c/d=(ac)/(bd)
    """
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    K=parent(frac1[0])
    if frac1[0]==0:
        return [K(0),K(1)]
    elif frac2[0]==0:
        return [K(0),K(1)]
    elif frac1[0]==frac2[0] and frac1[1]==frac2[1]:
        return [frac1[0]**2,frac1[1]**2]
    elif frac1[0]==frac2[0]:
        return [frac1[0]**2,frac1[1]*frac2[1]]
    elif frac1[1]==frac2[1]:
        return [frac1[0]*frac2[0],frac1[1]**2]
    else:
        return [frac1[0]*frac2[0],frac1[1]*frac2[1]]



def Frac_div(frac1:list,frac2:list):
    """ 
    Division of fractions.
    For given two fractions a/b, c/d, 
    give (a/b)/(c/d)=(ad)/(bc)
    """
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    return [frac1[0]*frac2[1],frac1[1]*frac2[0]]




def Projective_Theta(frac_coord:list):
    """ 
    From theta coordinate with different denominators, 
    output "projective" theta coordinate as 4 elements.
    """
    assert(len(frac_coord)==4)
    [[n1,d1],[n2,d2],[n3,d3],[n4,d4]]=frac_coord
    d12=d1*d2
    d34=d3*d4
    return [n1*d2*d34,d1*n2*d34,d12*n3*d4,d12*d3*n4]




def Multpower_straight(value,max_exp:int):
    """ 
    For x, and integer N, 
    compute [1,x,x^2,x^3,...,x^N].
    """
    pow_list=[parent(value)(1),value]
    value_pow=value
    for exp in range(2,max_exp+1):
        value_pow*=value
        pow_list.append(value_pow)   
    assert(len(pow_list)==max_exp+1) 
    return pow_list





def Multpower_sq(value,l:int):
    """ 
    For x, and integer N, 
    compute [1,x,x^4,x^9,...,x^(N^2)].
    """
    max_index_index=2*len(ZZ(l-1).digits(2))-1
    value_2_pow=[value]
    for index_index in range(1,max_index_index+1):
        a=value_2_pow[index_index-1]
        value_2_pow.append(a**2)
    reuslt_dict={0**2:1,1**2:value}
    for k in range(2,l):
        bin_exp_ksq=ZZ(k**2).digits(2)
        pow_value=1
        for i in range(0,len(bin_exp_ksq)):
            if bin_exp_ksq[i]==1:
                pow_value*=value_2_pow[i]
        reuslt_dict[k**2]=pow_value
    return reuslt_dict
            
     
    
    
    