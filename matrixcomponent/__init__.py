JSON_VERSION = 7  # This version added adjacent connectors in departures column

import logging

LOGGER = logging.getLogger(__name__)

LOGGING_FORMAT_STR = \
    '[%(asctime)s - %(levelname)s - %(name)s - %(lineno)d] %(message)s'
"""str: format string for log messages"""

LOGGING_DATE_FORMAT = '%d/%m/%Y %H:%M:%S'
"""str: format string for time/date stamps in log messages."""

if len(LOGGER.handlers) == 0:
    logging.basicConfig(level=logging.INFO,
                        format=LOGGING_FORMAT_STR,
                        datefmt=LOGGING_DATE_FORMAT)
