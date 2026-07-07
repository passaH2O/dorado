import pytest
import logging
import warnings
from dorado import logging_config

def test_setup_logging():
    # Test default
    logging_config.setup_logging()
    assert logging_config.logger.level == logging.WARNING

    # Test setting level
    logging_config.setup_logging("INFO")
    assert logging_config.logger.level == logging.INFO

    # Test invalid level
    with pytest.raises(ValueError, match="Invalid logging level"):
        logging_config.setup_logging("INVALID")

def test_handle_verbose_deprecation():
    # Test verbose=True
    with pytest.warns(DeprecationWarning, match="The 'verbose' parameter is deprecated"):
        logging_config.handle_verbose_deprecation(True)
    assert logging_config.logger.level == logging.INFO

    # Test verbose=False
    with pytest.warns(DeprecationWarning, match="The 'verbose' parameter is deprecated"):
        logging_config.handle_verbose_deprecation(False)
    assert logging_config.logger.level == logging.ERROR

    # Test verbose=None
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        logging_config.handle_verbose_deprecation(None)
        # Assuming no other warnings are triggered
        assert not any(issubclass(warn.category, DeprecationWarning) for warn in w)
