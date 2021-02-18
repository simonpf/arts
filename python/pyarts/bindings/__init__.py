"""
pyarts.bindings
===============

The pyarts.bindings module contains provides access to native ARTS classes
defined in C++ via bindings generated using pybind11. The top-level module
also provides some convenience functions to determine whether a given
object corresponds to a class exposed via these bindings.
"""
import re
import sys

from pyarts.bindings import cxx

_RE_MATCH_LOW_TO_CAP = re.compile(r"(\w)([A-Z][a-z0-9]+)")
_RE_MATCH_CAP_TO_LOW = re.compile(r"([a-z0-9]+)([A-Z])")

def _to_snake(camel):
    """
    Turns name given in CamelCase into snake_case.

    Args:
        The input name in CamelCase.

    Return:
        The given name transcribed to snake_case.
    """
    camel = re.sub(_RE_MATCH_LOW_TO_CAP, r"\1_\2", camel)
    return re.sub(_RE_MATCH_CAP_TO_LOW, r"\1_\2", camel).lower()

def has_bindings(value):
    """
    Determines whether a given value is a Python representation of a native
    ARTS class by checking whether parent module provides a copy_to_wsv method.
    """
    if hasattr(value, "__module__"):
        value_module = sys.modules[value.__module__]
        return hasattr(value_module, "copy_to_wsv")
    return None
    return hasattr(value_module, "copy_to_wsv")

def get_bindings_module(value):
    """
    Return:

    The sub-module containing the bindings for the class corresponding
    to the given values or None if value is not an instance of a class
    exposed through the C++ bindings.
    """
    return sys.modules.get(value.__module__)

def get_bindings_module_for_group(group):
    """
    Return:

    The sub-module containing the bindings for the given WSV group
    None if the corresponding class is not exposed through the C++
    bindings.
    """
    return getattr(cxx, _to_snake(group), None)

def get_class_for_group(group):
    """
    Return the class providing bindings for the given group.

    Args:
        group: The name of the ARTS group for which to return the
            corresponding class.

    Return:
        class object implementing the bindings for the given group
        or None.
    """
    module = get_bindings_module_for_group(group)
    if not module:
        return None
    return getattr(module, group)
