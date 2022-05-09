import time


def decoratortimer(decimal, logger):
    """get the execution time of a function

    Args:
        decimal (int): decimal precision of evaluated time
    """
    def decoratorfunction(f):
        def wrap(*args, **kwargs):
            time1 = time.monotonic()
            result = f(*args, **kwargs)
            time2 = time.monotonic()
            msg = "function '{:s}' took {:.{}f} s".format(
                f.__name__, ((time2-time1)), decimal)
            logger.info(msg)
            return result
        return wrap
    return decoratorfunction


def mean(a):
    """calculate mean

    Args:
        a (list): input list of values

    Returns:
        float: mean of the list
    """
    return sum(a)/len(a)
