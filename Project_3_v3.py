import math
import matplotlib.pyplot as plt
import scipy.integrate as intg
import numpy as np


# Global Variables
rho_max = 1.
u_max = 1.
x_center = .5
lam = .1
n = 150
perturbation = 10.**(-3.)


def init_grid(delta_rho,rho_bar):

    h = 1./n
    
    cells = []
    xs = []

    for i in range(n):
        x = (i*h)+(h/2.)
        exp_term = math.exp((-(x-x_center)**2.)/(lam**2.))

        if (i == n-1):
            cells.append(cells[0])
        else:
            cells.append(rho_bar+(delta_rho*exp_term))

        xs.append(x)

    #plt.plot(xs,cells)
    plt.show()

    return cells, xs


def f(rho):
    return (rho*u_max)*(1.-(rho/rho_max))


def rho_system(cells,t):
    
    #print(cells)

    n = len(cells)
    h = 1./n

    df = []
    drho = []

    for i in range(n):
        if (i == n-1):
            f_plus = f(cells[0])
            f_minus = f(cells[i-1])

            dfdx = (f_plus-f_minus)/(2.*h)
            drhodt = -dfdx

            df.append(dfdx)
            drho.append(drhodt)

        else:
            f_plus = f(cells[i+1])
            f_minus = f(cells[i-1])

            dfdx = (f_plus-f_minus)/(2.*h)
            drhodt = -dfdx

            df.append(dfdx)
            drho.append(drhodt)
            
    return drho



def SHOCK_CAPTURE(cells, t):
    drho = []

    
    #Reconstruction step
def minmod(a,b):
    
    if (a*b > 0):
        if abs(a)>abs(b):
            val = b
        else:
            val = a
    else:
        val = 0
    
    return val

    
def U_Left(Ui, Ui_plus_1, Ui_minus_1):
    a = Ui_plus_1 - Ui
    b = Ui - Ui_minus_1
    
    UL = Ui + 0.5*minmod(a,b)
    
    return UL


def U_Right(Ui, Ui_plus_1, Ui_plus_2):
    a = Ui_plus_2 - Ui_plus_1
    b = Ui_plus_1 - Ui
    
    UR = Ui_plus_1 - 0.5*minmod(a,b)
    
    return UR

    #Solving Riemann Probem step

    
def check_S(rho_L, rho_R):
    
    f_of_rho_L = f(rho_L)
    #print("rho_L", rho_L)
    #print("f_of_rho_L", f_of_rho_L)
    f_of_rho_R = f(rho_R)
    #print("rho_R", rho_R)
    #print("f_of_rho_R", f_of_rho_R)
    
    S = (f_of_rho_L - f_of_rho_R) / (rho_L - rho_R)
    #print("S", S)
    
    if S > 0:
        fi_plus_one_half = f_of_rho_L
    elif S < 0:
        fi_plus_one_half = f_of_rho_R
    else:
        fi_plus_one_half = f_of_rho_R                  #this is for the n-1 term, where S will equal 0/0
    
    return fi_plus_one_half


def METH(cells, t):
    
    #print("cells", cells)
    #print("t", t)
    n = len(cells)
    #print("length", n)
    h = 1./n
    f_i = []
    drho = []
    

    for i in range(n):
        if (i == n-1):
            U__i = cells[i]
            U__i_minus_one = cells[i-1]
            U__i_plus_one = cells[0]
            U__i_plus_two = cells[1]
            
            #print("Ui for n-1", U__i)
            #print("U i-1 for n-1", U__i_minus_one)
            #print("U i+1 for n-1", U__i_plus_one)
            #print("U i+2 for n-1", U__i_plus_two)
            
            
        elif (i == n-2):
            U__i = cells[i]
            U__i_minus_one = cells[i-1]
            U__i_plus_one = cells[i+1]
            U__i_plus_two = cells[0]
            
            #print("Ui for n-2", U__i)
            #print("U i-1 for n-2", U__i_minus_one)
            #print("U i+1 for n-2", U__i_plus_one)
            #print("U i+2 for n-2", U__i_plus_two)
            
        else:
            U__i = cells[i]
            U__i_minus_one = cells[i-1]
            U__i_plus_one = cells[i+1]
            U__i_plus_two = cells[i+2]
            
            #print("Ui for", i, U__i)
            #print("U i-1 for", i, U__i_minus_one)
            #print("U i+1 for", i, U__i_plus_one)
            #print("U i+2 for", i, U__i_plus_two)
        
            
        ULeft = U_Left(U__i, U__i_plus_one, U__i_minus_one)
        #print("UL", ULeft)
        URight = U_Right(U__i, U__i_plus_one, U__i_plus_two)
        #print("UR", URight)
        
        f_i.append(check_S(ULeft, URight))
        
    for j in range(n):
        
        if (j == n-1):
            f_iplus = f_i[0]
            f_iminus = f_i[j-1]

            dfdx = (f_iplus-f_iminus)/(h)
            drhodt = -dfdx
            
            drho.append(drhodt)

        else:
            f_iplus = f_i[j+1]
            f_iminus = f_i[j-1]

            dfdx = (f_iplus-f_iminus)/(h)
            drhodt = -dfdx

            drho.append(drhodt)
    
    return drho
        
        
    
def integrate_d(func,rhos):
    t = range(n)
    solution = intg.odeint(func,rhos,t)
    return solution

def rho_plots():
    cells, x_vals = init_grid(rho_max*(10**-2),rho_max*.6)
    rho_t = integrate_d(METH,cells)

    for num in range(5):
        maximum = max(rho_t[num])
        print("y-max at t", num, ":", maximum)
        counter = 0
        for val in rho_t[num]:
            if val == maximum:
                x_max = x_vals[counter]
            else:
                counter += 1
        print("x-max at t", num, ":", x_max, "\n")

    plt.plot(x_vals,rho_t[0],label='t0')
    plt.plot(x_vals,rho_t[5],label='t5')
    plt.plot(x_vals,rho_t[10],label='t10')
    plt.plot(x_vals,rho_t[20],label='t20')
    plt.plot(x_vals,rho_t[45],label='t45')
    plt.plot(x_vals,rho_t[60],label='t60')
    plt.plot(x_vals,rho_t[75],label='t75')
    #plt.plot(x_vals,rho_t[100],label='t100')
    plt.xlabel('x')
    plt.ylabel('rho')
    plt.title('Rho vs. position (with rho_bar=rho_max/2, perturbation=.001*rho_max)')
    plt.legend(loc='best')
    plt.show()

rho_plots()



#integrate_d(rho_system,init_grid(10.**-3.,rho_max/2.)[0])
