import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import itertools

afivo = np.loadtxt('Afivo_data.dat', delimiter=',', skiprows=0)
zdplaskin = np.loadtxt('ZDPlaskin_simplerAir.dat', delimiter=' ', skiprows=0)


a_species = 'e electric_fld photo N2_plus O O2_min O2_plus O_min'.split()
a_species.remove('electric_fld')
a_species.remove('photo')
afivo = np.delete(afivo, 2, 1)
afivo = np.delete(afivo, 2, 1)
colors = iter(cm.rainbow(np.linspace(0, 1, len(a_species))))
for meh in range(len(a_species)):
  #plt.yscale('log')
  blah = next(colors)
  plt.plot(afivo[:,0], afivo[:,meh+1], label=a_species[meh], color=blah, marker='.')
  #
  plt.plot(zdplaskin[:,0],1e6*zdplaskin[:,meh+1], label=a_species[meh], color = blah, marker='v')
  plt.yscale('log')
  plt.legend()
#plt.plot(zdplaskin[:,0], zdplaskin[:,2])
plt.show()
