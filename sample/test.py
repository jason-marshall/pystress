from pystress import *
import matplotlib.pyplot as plt

# main() function
def main():

    s11 = 8.0
    s12 = -2.0
    s22 = 6.0

    stress = StressState(s11=s11,s22=s22,s12=s12,symmetric=Symmetric.yes)
    stress.plot_stress()
    stress.plot_mohr_coulomb(c=0.0, phi=20.0)

    fig = plt.gcf()
    fig.set_size_inches(7.0, 5.0)

    plt.savefig('test-7.eps', dpi=600, format='eps')

    plt.show()
    
# call main
if __name__ == '__main__':
  main()
