# -*- coding: utf-8 -*-
"""
Logging configuration for the dorado package.

Provides centralized logging setup and backward-compatible handling
of the deprecated `verbose` parameter.

Project Homepage: https://github.com/passaH2O/dorado
"""
import logging
import warnings

# Package-level logger
logger = logging.getLogger("dorado")


def setup_logging(level="WARNING", handler=None):
    """Configure dorado logging.

    **Inputs** :

        level : `str`, optional
            Logging level. One of "DEBUG", "INFO", "WARNING", "ERROR",
            "CRITICAL". Default is "WARNING".

        handler : `logging.Handler`, optional
            Custom logging handler. If not provided, a StreamHandler
            writing to stderr is used.

    **Example**::

        import dorado

        # Show informational messages
        dorado.setup_logging("INFO")

        # Log to a file
        import logging
        dorado.setup_logging("DEBUG", handler=logging.FileHandler("dorado.log"))

    """
    logger.setLevel(getattr(logging, level.upper()))
    if not logger.handlers:
        h = handler or logging.StreamHandler()
        h.setFormatter(
            logging.Formatter("%(name)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(h)


def handle_verbose_deprecation(verbose):
    """Map deprecated verbose parameter to logging level.

    This function provides backward compatibility for the deprecated
    `verbose` parameter. It configures the logger based on the verbose
    value and emits a deprecation warning.

    **Inputs** :

        verbose : `bool` or `None`
            Legacy verbose parameter. If True, sets logging to INFO level.
            If False, sets logging to ERROR level. If None, does nothing.

    """
    if verbose is not None:
        warnings.warn(
            "The 'verbose' parameter is deprecated. "
            "Use dorado.setup_logging() instead.",
            DeprecationWarning,
            stacklevel=3
        )
        level = "INFO" if verbose else "ERROR"
        setup_logging(level)
