import sys
import inspect
from numba import jit
from functools import wraps

def old_sphinx_compat_jit(f_or_has_args=None, **kwargs):
    """To be used instead of @jit in order to prevent jit from rewriting
    obscure parts of the function's documentation (to contain e.g. unquoted
    asterisks) when sphinx is trying to read it.

    should not change jit functionlity if you are not sphinx."""
    def jit_passthru(f):
        return jit(f, **kwargs)
    def do_nothing(f):
        return f
    # returns one "FrameInfo" object per stack frame
    stack = inspect.stack()
    # get the "name" of each frame's calling module, as in
    # https://stackoverflow.com/questions/7871319/how-to-know-who-is-importing-me-in-python?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    modules = [finfo.frame.f_globals.get('__name__') for finfo in stack]
    # if sphinx_compat_jit was used as a decorator with no parameters
    # or as a regular function
    if f_or_has_args is not None:
        if 'sphinx' in modules:
            return f_or_has_args
        else:
            return jit(f_or_has_args, **kwargs)
    # otherwise it was used as a decorator with parameters
    else:
        if 'sphinx' in modules:
            actual_decorator = do_nothing
        else:
            actual_decorator = jit_passthru
    return actual_decorator


def sphinx_compat_jit(f=None, **kwargs):
    """Fixes jit to correctly @wraps the function.

    Can be used as a drop-in replacement for jit.

    We want both of the following to work:

    >>> new_f_maker = sphinx_compat_jit(**kwargs)
    >>> f = new_f_maker(f)
    >>> # and
    >>> f = sphinx_compat_jit(f)

    Which correspond to

    >>> @sphinx_compat_jit(**kwargs)
    >>> def f():
    >>>     pass
    >>> # and (resp)
    >>> @sphinx_compat_jit
    >>> def f():
    >>>     pass
    """
    def jit_passthru(f):
        # if we were passed kwargs (like @jit(nopython=True))
        if f is None or not callable(f):
            jitter = jit(**kwargs)
        else:
            jitter = jit
        # now return the new f, which should just be jitter(f)!
        jitted_f = jitter(f)
        @wraps(f)
        def new_f(*args, **kw):
            return jitted_f(*args, **kw)
        return new_f
    # if we're being called as a decorator with no kwargs, like
    # @sphinx_compat_jit
    # def f():...
    # this is syntax sugar for
    # f = sphinx_compat_jit(f)
    # first detect that this is the case...
    if f is not None and callable(f):
        # now return "new_f"
        return jit_passthru(f)
    # if we're being called with kwargs, like
    # @sphinx_compat_jit(nopython=True)
    # def f():...
    # this is syntax sugar for
    # jitter = jit(no_python=True)
    # f = jitter(f)
    else:
        # so return the function that "makes" new_f
        return jit_passthru

