import FiberBundleModel_constants as FBM_c
#import FiberBundleModel_ModelParameters as FBM_MP
import numpy as np
import matplotlib.pyplot as plt





def matrix_k(b_k):
    m =[[6,0,0,b_k],[0,2,0,0],[0,0,1,0],[0,0,0,1]]
    return m

def matrix_kMinus1():
    m = [[-6,0,0,0],[-6,-2,0,0],[-3,-2,-1,0],[-1,-1,-1,-1]]
    return m


def makeAMatrix(systemSize, b_k):
    A = np.zeros((4*systemSize, 4*systemSize))
    for k in range(systemSize):
        A[4*k:4*(k+1),4*k:4*(k+1)] = matrix_k(b_k) # Usikker på om matrisen skrives som vi vil
        if (k==0):
            A[0:4,3*systemSize:4*systemSize] = matrix_kMinus1()
        else:
            A[4*k:4*(k+1),4*(k-1):4*k] = matrix_kMinus1()
    return A

def makeuMatrix(systemSize, alfa):
    u = np.zeros(systemSize * 4)
    for k in range(systemSize):
        u[4*k:4*(k+1)] = [-alfa, -alfa / 2, -alfa / 6, -alfa / 24]
    return u
def solve(N, alfa, b_k):
    u = makeuMatrix(N, alfa)
    A = makeAMatrix(N, b_k)
    #print(A)
    #print(u)
    x = np.linalg.solve(A, u)
    #print(x)
    return x
#solve(4, FBM_c.alpha(9.81,1,0.5,FBM_c.B(0.1,2*10**11)), FBM_c.beta(1,0.01,0.1,0.5,FBM_c.B(0.1,2*10**11),4))

#PLOTTING



def plotSolution(solutionVector,numberOfSolutionIntervals, alpha):

    fig_solution=plt.figure(1,[16,10])

    fig2=plt.figure(2,[16,10])

    fig_solution.suptitle('$\eta$ og dens deriverte',fontsize=25)

    ax_y=fig2

    #ax_y=fig_solution.add_subplot(221)
    ax_dy = fig_solution.add_subplot(222)
    ax_ddy = fig_solution.add_subplot(223)
    ax_dddy = fig_solution.add_subplot(224)
    # fig1=plt.figure(1, [10, 5])
    # ax_y=fig1.add_subplot(221)
    #
    # fig2= plt.figure(2, [10, 5])
    # ax_dy=fig2.add_subplot(221)
    #
    # fig3 = plt.figure(3, [10, 5])
    # ax_ddy=fig3.add_subplot(221)
    #
    # fig4 = plt.figure(4, [10, 5])
    # ax_dddy = fig4.add_subplot(221)

    fig2.suptitle("$\eta$")
    #ax_y.set_title("$\eta(ksi)$",fontsize=20)
    ax_dy.set_title("$d\eta/dksi$",fontsize=20)
    ax_ddy.set_title("$d^2\eta/dksi^2$",fontsize=20)
    ax_dddy.set_title("$d^3\eta/dksi^3$",fontsize=20)

    fig2.set_xlim([0,1])
    #ax_y.set_xlim()
    ax_dy.set_xlim([0, 1])
    ax_ddy.set_xlim([0, 1])
    ax_dddy.set_xlim([0, 1])

    # ax_y.set_ylabel("$\eta$")
    # ax_dy.set_ylabel("$d\eta/dksi$")
    # ax_ddy.set_ylabel("$dˆ2\eta/dksiˆ2$")
    # ax_dddy.set_ylabel("$dˆ3\eta/dksiˆ3$")

    xi_minus_k = np.linspace(0.0, 1, 200)
    for n in range(len(alpha)):
        alph=alpha[n]
        solVector=solutionVector[n]
        for k in range(numberOfSolutionIntervals): # Dette blir vel riktig


            eta = -(alph / 24) * xi_minus_k ** 4 + solVector[4*k] * xi_minus_k ** 3 + solVector[4*k+1] * xi_minus_k ** 2 + solVector[4*k+2]*xi_minus_k+solVector[4*k+3]



            d_eta=-(alph/6)*xi_minus_k**3+3*solVector[4*k]*xi_minus_k**2+2*solVector[4*k+1]*xi_minus_k+solVector[4*k+2]

            dd_eta=-(alph/2)*xi_minus_k**2 + 6*solVector[4*k]*xi_minus_k + 2*solVector[4*k+1]

            ddd_eta= -alph*xi_minus_k + 6*solVector[4*k]

            fig2.plot(xi_minus_k,eta)
            #ax_y.plot(xi_minus_k, eta)
            ax_dy.plot(xi_minus_k, d_eta)
            ax_ddy.plot(xi_minus_k,dd_eta)
            ax_dddy.plot(xi_minus_k,ddd_eta)

        #ax_y.legend(loc='upper center',labels="5")

    plt.show(fig_solution)
    plt.show(fig2)

plotSolution([solve(4,5.29e-09,53.05)], 4, [5.29e-09])

#Løsning av likning (18) for DNA-molekyl


# A_dna=makeAMatrix(4,53.05)
# u_dna=makeuMatrix(4,5.29e-09)
# x_dna=np.linalg.solve(A_dna,u_dna)
# print (x_dna)