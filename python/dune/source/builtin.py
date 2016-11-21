from __future__ import division, print_function, unicode_literals


class BuiltInFunction:
    def __init__(self, header, cppType, name, namespace='std', targs=None, args=None):
        self.header = header
        self.cppType = cppType
        self.name = name
        self.namespace = namespace
        self.tarts = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() for a in args]

    def __call__(self, *args):
        if len(args) != len(self.args):
            raise Exception('Wrong number of Arguments: ' + len(args) + ' (should be ' + len(self.args) + ').')
        from .expression import Application, makeExpression
        return Application(self, args=[makeExpression(arg) for arg in args])


atan = BuiltInFunction('cmath', 'X', 'atan', targs=['class X'], args=['conat X &x'])
atan2 = BuiltInFunction('cmath', 'X', 'atan2', targs=['class X'], args=['const X &x', 'const X &y'])
cos = BuiltInFunction('cmath', 'X', 'cos', targs=['class X'], args=['const X &x'])
pow_ = BuiltInFunction('cmath', 'X', 'pow', targs=['class X'], args=['const X &x', 'const X &y'])
sin = BuiltInFunction('cmath', 'X', 'sin', targs=['class X'], args=['const X &x'])
tan = BuiltInFunction('cmath', 'X', 'tan', targs=['class X'], args=['const X &x'])

def get(i):
    return BuiltInFunction('tuple', 'auto', 'get< ' + str(i) + ' >', targs=['class T'], args=['const T &arg'])

make_pair = BuiltInFunction('utility', 'std::pair< U, V >', 'make_pair', targs=['class U', 'class V'], args=['const U &left', 'const V &right'])