import time


def decoratortimer(decimal):
    """get the execution time of a function

    Args:
        decimal (int): decimal precision of evaluated time
    """
    def decoratorfunction(f):
        def wrap(*args, **kwargs):
            time1 = time.monotonic()
            result = f(*args, **kwargs)
            time2 = time.monotonic()
            print('{:s} function took {:.{}f} ms'.format(
                f.__name__, ((time2-time1)*1000.0), decimal))
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