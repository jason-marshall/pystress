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
    fig.set_size_inches(3.4, 2.5)

    plt.savefig('test-double-column.ps', dpi=600, format='ps')

    plt.show()
    
# call main
if __name__ == '__main__':
  main()
