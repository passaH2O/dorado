"""Plot for quickstart demo 1."""
import matplotlib.pyplot as plt
from particlerouting import routines
from particlerouting.example_data.define_params import make_rcm_params

params = make_rcm_params()
plt.imshow(params.depth)
plt.title('Water Depth')
plt.show()
