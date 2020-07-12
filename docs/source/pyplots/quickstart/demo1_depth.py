"""Plot for quickstart demo 1."""
import matplotlib.pyplot as plt
from particlerouting.example_data.define_params import make_rcm_params

params = make_rcm_params()
plt.imshow(params.depth)
plt.colorbar(fraction=0.025, pad=0.03)
plt.title('Water Depth')
plt.show()
