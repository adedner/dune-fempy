from __future__ import division, print_function, unicode_literals

from .common import Block, Statement
from .operator import BinaryOperator, BracketOperator, PrefixUnaryOperator, PostfixUnaryOperator


class Expression(Statement):
    def __init__(self, cppType=None):
        Statement.__init__(self)
        self.cppType = cppType

    def __add__(self, other):
        return Application(BinaryOperator('+'), args=(self, makeExpression(other)))

    def __sub__(self, other):
        return Application(BinaryOperator('-'), args=(self, makeExpression(other)))

    def __mul__(self, other):
        return Application(BinaryOperator('*'), args=(self, makeExpression(other)))

    def __truediv__(self, other):
        return Application(BinaryOperator('/'), args=(self, makeExpression(other)))

    def __mod__(self, other):
        return Application(BinaryOperator('%'), args=(self, makeExpression(other)))

    def __iadd__(self, other):
        return Application(BinaryOperator('+='), args=(self, makeExpression(other)))

    def __isub__(self, other):
        return Application(BinaryOperator('-='), args=(self, makeExpression(other)))

    def __imul__(self, other):
        return Application(BinaryOperator('*='), args=(self, makeExpression(other)))

    def __imod__(self, other):
        return Application(BinaryOperator('%='), args=(self, makeExpression(other)))

    def __itruediv__(self, other):
        return Application(BinaryOperator('/='), args=(self, makeExpression(other)))

    def __lt__(self, other):
        return Application(BinaryOperator('<'), args=(self, makeExpression(other)))

    def __le__(self, other):
        return Application(BinaryOperator('<='), args=(self, makeExpression(other)))

    def __eq__(self, other):
        return Application(BinaryOperator('=='), args=(self, makeExpression(other)))

    def __ne__(self, other):
        return Application(BinaryOperator('!='), args=(self, makeExpression(other)))

    def __ge__(self, other):
        return Application(BinaryOperator('>='), args=(self, makeExpression(other)))

    def __gt__(self, other):
        return Application(BinaryOperator('>'), args=(self, makeExpression(other)))

    def __neg__(self):
        return Application(PrefixUnaryOperator('-'), args=(self,))

    def __pos__(self):
        return Application(PrefixUnaryOperator('+'), args=(self,))

    def __getitem__(self, index):
        if self.cppType is not None and self.cppType.startswith('std::tuple'):
            from .builtin import get
            return Application(get(index), args=(self,))
        else:
            if isinstance(index, tuple):
                result = self
                for i in index:
                    result = Application(BracketOperator(), args=(result, makeExpression(i)))
                return result
            else:
                return Application(BracketOperator(), args=(self, makeExpression(index)))


class Application(Expression):
    def __init__(self, function, args=None):
        Expression.__init__(self)
        self.function = function
        self.args = tuple(args)

    def __hash__(self):
        return hash((self.cppType, self.function, self.args))


class ConstantExpression(Expression):
    def __init__(self, cppType, value):
        Expression.__init__(self, cppType)
        self.value = value

    def __hash__(self):
        return hash((self.cppType, self.value))


class ConstructExpression(Expression):
    def __init__(self, cppType, args=None):
        Expression.__init__(self, cppType)
        self.args = None if args is None else [makeExpression(arg) for arg in args]

    def __hash__(self):
        return hash((self.cppType, self.args))


class InitializerList(Expression):
    def __init__(self, *args):
        Expression.__init__(self)
        self.args = tuple(args)


class LambdaExpression(Expression):
    def __init__(self, args=None, capture=None, code=None):
        Expression.__init__(self, None)
        self.args = None if args is None else tuple(args)
        self.capture = capture
        if code is None:
            self.code = None
        elif isinstance(code, Block):
            self.code = tuple(block.content)
        elif isinstance(code, (list, set, tuple)):
            self.code = (o for o in code)
        else:
            self.code = (code,)

    def __hash__(self):
        return hash((self.cppType, self.args, self.capture, self.code))


class Variable(Expression):
  def __init__(self, cppType, name):
      Expression.__init__(self, cppType)
      self.name = name

  def __hash__(self):
      return hash((self.cppType, self.name))


def makeExpression(expr):
    if isinstance(expr, bool):
        return ConstantExpression('bool', 'true' if expr else 'false')
    elif isinstance(expr, int):
        if expr < 0:
            return -ConstantExpression('int', str(-expr))
        else:
            return ConstantExpression('int', str(expr))
    elif isinstance(expr, float):
        s = str(abs(expr))
        if "." not in s:
            s += ".0"
        e = ConstantExpression('double', s)
        return -e if expr < 0 else e
    else:
        return expr


def assign(left, right):
    return Application(BinaryOperator('='), args=(left, makeExpression(right)))


def construct(cppType, *args):
    return ConstructExpression(cppType, args if args else None)


def lambda_(args=None, capture=None, code=None):
    return LambdaExpression(args=args, capture=capture, code=code)