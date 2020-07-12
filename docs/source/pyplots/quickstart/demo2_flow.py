"""Plot for quickstart demo 2."""
import matplotlib.pyplot as plt
from particlerouting.example_data.define_params import make_anuga_params

anugaparams = make_anuga_params()
plt.figure(figsize=(8, 5))
plt.subplot(1, 2, 1)
plt.imshow(anugaparams.qx)
plt.colorbar(fraction=0.05, pad=0.03)
plt.title('x-component of water discharge')
plt.subplot(1, 2, 2)
plt.imshow(anugaparams.qy)
plt.colorbar(fraction=0.05, pad=0.03)
plt.title('y-component of water discharge')
plt.tight_layout()
plt.show()
