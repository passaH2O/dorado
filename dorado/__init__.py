from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pydorado")
except PackageNotFoundError:
    __version__ = "unknown"


from . import lagrangian_walker
from . import parallel_routing
from . import particle_track
from . import routines
from . import spatial
from .logging_config import setup_logging, logger
