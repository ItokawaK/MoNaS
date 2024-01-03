import logging
from logging import getLogger, StreamHandler, FileHandler, Formatter

logger = getLogger(__name__)

def set(log_file):
    logger.parent.setLevel(logging.DEBUG)
    stream_handler = StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    handler_format = Formatter('[%(levelname)s: %(asctime)s] %(message)s (%(name)s)')
    stream_handler.setFormatter(handler_format)

    file_handler = FileHandler(log_file, 'a')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(handler_format)
    logger.parent.addHandler(stream_handler)
    logger.parent.addHandler(file_handler)
