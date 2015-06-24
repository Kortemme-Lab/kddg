#!/usr/bin/python2.4
# encoding: utf-8
"""
api_layers.py
The definition of the layers of the database API and the generic user interface class.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""


import inspect
import functools
from tools import colortext

### API function decorators. These are used to group functions together when printing the help text.

functional_layer = {
    0 : 'API warnings',
    1 : 'Information layer',
    2 : 'Prediction layer',
    3 : 'Results layer',
    4 : 'Analysis layer',
    5 : 'Application layer',
    None: 'Miscellanous'
}

def alien(func):
    func._helptype = 'Alien functions (these should be moved into another package)'
    func._layer = 0
    func._layer_order = 0
    return func

def brokenfn(func):
    func._helptype = 'Broken functions: this need to be fixed/updated'
    func._layer = 0
    func._layer_order = 1
    return func

def deprecated(func):
    func._helptype = 'Deprecated functions'
    func._layer = 0
    func._layer_order = 2
    return func

def informational_misc(func):
    func._helptype = 'Miscellaneous information API'
    func._layer = 1
    func._layer_order = 0
    return func

def informational_file(func):
    func._helptype = 'File information API'
    func._layer = 1
    func._layer_order = 1
    return func

def informational_pdb(func):
    func._helptype = 'Structure information API'
    func._layer = 1
    func._layer_order = 2
    return func

def informational_job(func):
    func._helptype = 'Prediction information API'
    func._layer = 1
    func._layer_order = 3
    return func

def job_creator(func):
    func._helptype = 'Job creation API'
    func._layer = 2
    func._layer_order = 0
    return func

def job_input(func):
    func._helptype = 'Input file generation API'
    func._layer = 2
    func._layer_order = 1
    return func

def job_execution(func):
    func._helptype = 'Job execution API'
    func._layer = 2
    func._layer_order = 2
    return func

def job_completion(func):
    func._helptype = 'Job completion API'
    func._layer = 2
    func._layer_order = 3
    return func

def job_results(func):
    func._helptype = 'Results API'
    func._layer = 3
    func._layer_order = 0
    return func

def analysis_api(func):
    func._helptype = 'Analysis API'
    func._layer = 4
    func._layer_order = 0
    return func

def app_pymol(func):
    func._helptype = 'PyMOL API'
    func._layer = 5
    func._layer_order = 0
    return func


class GenericUserInterface(object):
    '''This is the class that should be used to interface with the database. It hides functions that should only be called
       within this other API functions.

       The class contains a private copy of the internal API and wraps the public functions of that API so that the
       functions of GenericUserInterface contain only the public functions of the internal API. Private functions
       are denoted as such by a leading underscore in the function name.
       '''

    @staticmethod
    def generate(cls, passwd = None, username = 'kortemmelab'):
        return GenericUserInterface(cls, passwd = passwd, username = username)

    @staticmethod
    def bind_object_function(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs): return fn(*args, **kwargs)
        return wrapper

    def __init__(self, cls, passwd = None, username = 'kortemmelab'):

        self._ddg_interface = cls(passwd = passwd, username = username)
        self._api_functions = []
        self._api_function_args = {}
        self.DDG_db = self._ddg_interface.DDG_db
        self.DDG_db_utf = self._ddg_interface.DDG_db_utf
        self.cls = cls

        for m in inspect.getmembers(cls, predicate=inspect.ismethod):
            if m[0][0] != '_':
                fn_name = m[0]
                fn_ref = getattr(self._ddg_interface, fn_name)
                self._api_function_args[fn_name] = fn_ref.func_code.co_varnames[:fn_ref.func_code.co_argcount]
                self._api_functions.append(fn_name)
                self.__dict__[fn_name] = GenericUserInterface.bind_object_function(getattr(self._ddg_interface, fn_name))


    def help(self, show_deprecated_functions = False):
        print(self.get_help(show_deprecated_functions = show_deprecated_functions))


    def get_help_str(self, fn, fn_name):

        if fn.__doc__:
            return fn.__doc__

        # Wind up the hierarchy until we find the class where this function was last defined
        for parent in self.cls.__mro__[1:]:
            overridden = getattr(parent, fn_name, None)
            if overridden: break
        return overridden.__doc__


    def get_help(self, show_deprecated_functions = False):
        helpstr = []
        title = ' %s API ' % self._ddg_interface.__class__.__name__
        l = len(title)
        helpstr.append(colortext.mcyan('\n' + ('*' * (l + 10)) + '\n' + ('*' * 5) + title + ('*' * 5) + '\n' + ('*' * (l + 10)) + '\n'))

        doc_strings = {}
        for fn_name in sorted(self._api_functions):
            fn = self.__dict__[fn_name]

            function_layer, function_layer_order, function_class = None, None, None
            try:
                function_layer = fn._layer
                assert(function_layer in functional_layer)
                function_layer_order = fn._layer_order
            except:
                function_layer = None
                function_layer_order = 0
            try:
                function_class = fn._helptype
            except:
                function_class = 'Miscellanous'

            if function_class == 'Deprecated functions' and not show_deprecated_functions:
                continue

            doc_strings[function_layer] = doc_strings.get(function_layer, {})
            doc_strings[function_layer][function_layer_order] = doc_strings[function_layer].get(function_layer_order, {})
            doc_strings[function_layer][function_layer_order][function_class] = doc_strings[function_layer][function_layer_order].get(function_class, {})
            doc_strings[function_layer][function_layer_order][function_class][fn_name] = self.get_help_str(fn, fn_name)

        for function_layer, function_layer_components in sorted(doc_strings.iteritems()):
            function_layer_name = functional_layer[function_layer]
            prefix = ''
            if function_layer != None:
                prefix = 'Layer %d: ' % function_layer
            helpstr.append(colortext.mcyan('-------- %s%s --------\n' % (prefix, function_layer_name)))
            for function_layer_order, function_classes in sorted(function_layer_components.iteritems()):
                for function_class, fn_names in sorted(function_classes.iteritems()):
                    helpstr.append(colortext.mlightpurple('  %s\n' % function_class))
                    for fn_name, docstr in sorted(fn_names.iteritems()):
                        helpstr.append(colortext.mgreen('  %s(%s)' % (fn_name, ', '.join(self._api_function_args[fn_name]))))
                        if docstr:
                            helpstr.append(colortext.myellow('    %s' % ('\n    '.join([s.strip() for s in docstr.split('\n') if s.strip()]))))
                        else:
                            helpstr.append(colortext.mred('    <not documented>'))
                        helpstr.append('')
        return '\n'.join(helpstr)
