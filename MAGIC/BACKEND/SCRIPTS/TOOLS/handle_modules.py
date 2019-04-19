#https://stackoverflow.com/questions/1668223/how-to-de-import-a-python-module

import imp
import os
__MODULE_EXTENSIONS = ('.py', '.pyc', '.pyo')

def package_contents(package_name):
    file, pathname, description = imp.find_module(package_name)
    if file:
        raise ImportError('Not a package: %r', package_name)
    # Use a set because some may be both source and compiled.
    return set([os.path.splitext(module)[0]
        for module in os.listdir(pathname)
        if module.endswith(__MODULE_EXTENSIONS)])

def delete_package(packname):
    for _mod in package_contents(packname):
        if _mod != '__init__':
            delete_module(_mod, package=packname)

def delete_module(modname, package=None, paranoid=None):
    from sys import modules
    try:

        if package != None:
            thismod = modules[package + '.' + modname] # SM
        else:
            thismod = modules[modname]

    except KeyError:
        raise ValueError(modname)
    these_symbols = dir(thismod)
    if paranoid:
        try:
            paranoid[:]  # sequence support
        except:
            raise ValueError('must supply a finite list for paranoid')
        else:
            these_symbols = paranoid[:]

    if package != None:
        del modules[package + '.' + modname] # SM
    else:
        del modules[modname]

    for mod in modules.values():
        try:
            delattr(mod, modname)
        except AttributeError:
            pass
        if paranoid:
            for symbol in these_symbols:
                if symbol[:2] == '__':  # ignore special symbols
                    continue
                try:
                    delattr(mod, symbol)
                except AttributeError:
                    pass
