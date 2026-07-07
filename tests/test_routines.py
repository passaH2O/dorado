import pytest
import warnings
from dorado import routines

class DummyParams:
    def __init__(self):
        self.dx = 1
        self.depth = [[1]]
        self.qx = [[0]]
        self.qy = [[1]]
        self.theta = 1
        self.model = 'DeltaRCM'

def test_routines_verbose_deprecation():
    # Test steady_plots
    walk_data = {'x': [], 'y': []}
    with pytest.warns(DeprecationWarning, match="The 'verbose' parameter is deprecated"):
        try:
            routines.steady_plots(walk_data, 1, verbose=True)
        except Exception:
            pass # We just care about the warning
            
    with pytest.warns(DeprecationWarning, match="The 'verbose' parameter is deprecated"):
        try:
            routines.time_plots(walk_data, 1, verbose=True)
        except Exception:
            pass

