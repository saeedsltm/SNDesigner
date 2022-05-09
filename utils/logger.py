import logging


def myLogger(loggerFile, mode, level=logging.INFO):
    """Method to return a custom logger with the given name and level

    Args:
        loggerFile (_type_): log file output name
        mode (_type_): new file or append to an existing file
        level (_type_, optional): logging level. Defaults to logging.INFO.

    Returns:
        logger: a looger will be returned
    """
    logger = logging.getLogger(loggerFile)
    logger.setLevel(level)
    format_string = ("%(asctime)s %(message)s")
    log_format = logging.Formatter(format_string)
    file_handler = logging.FileHandler("%s.log" % (loggerFile), mode=mode)
    file_handler.setFormatter(log_format)
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.addHandler(file_handler)
    return logger
