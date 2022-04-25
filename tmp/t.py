import numpy as np
import proplot as pplt
##pplt.rc['font.family'] = "TeX Gyre Schola"
state = np.random.RandomState(51423)
data = 2 * (state.rand(100, 5) - 0.5).cumsum(axis=0)
fig = pplt.figure(suptitle='Single subplot')
ax = fig.subplot(xlabel='x axis', ylabel='y axis')
ax.plot(data, lw=2)

pplt.show()
