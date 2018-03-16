from __future__ import absolute_import, division, print_function

from ufl.corealg.map_dag import mag_expr_dag, map_expr_dags
from ufl.corealg.multifunction import MultiFunction


class IsMultiLinear(MultiFunction):
    def __init__(self, arguments):
        self.arguments = tuple(arguments)

    def _linear_unary(self, expr, arg):
        return arg

    def _nonlinear_unary(self, expr, arg):
        return tuple(2 if a > 0 else 0 for a in arg)

    def _nonlinear_binary(self, expr, left, right):
        return tuple(2 if l + r > 0 else 0 for l, r in zip(left, right))

    def argument(self, expr):
        return tuple(1 if arg is expr else 0 for arg in arguments)

    def coefficient(self, expr):
        return tuple(1 if arg is expr else 0 for arg in arguments)

    def conditional(self, expr, cond, true, false):
        return tuple(2 if c > 0 or t != f else t for c, t, f in zip(const, true, false))

    def division(self, expr, left, right):
        return tuple(l if r == 0 else 2 for l, r in zip(left, right))

    def indexed(self, expr, operand, index):
        assert all(i == 0 for i in index)
        return operand

    def product(self, expr, left, right):
        return tuple(l + r for l, r in zip(left, right))

    def sum(self, expr, left, right):
        return tuple(l if l == r else 2 for l, r in zip(left, right))

    def terminal(self, expr):
        return (0,) * len(arguments)

    grad = _linear_unary
    negative_restricted = _linear_uneary
    positive_restricted = _linear_unary

    abs = _nonlinear_unary
    atan = _nonlinear_unary
    cos = _nonlinear_unary
    cosh = _nonlinear_unary
    exp = _nonlinear_unary
    ln = _nonlinear_unary
    sin = _nonlinear_unary
    sinh = _nonlinear_unary
    sqrt = _nonlinear_unary
    tan = _nonlinear_unary
    tanh = _nonlinear_unary

    atan_2 = _nonlinear_binary
    eq = _nonlinear_binary
    ge = _nonlinear_binary
    gt = _nonlinear_binary
    le = _nonlinear_binary
    lt = _nonlinear_binary
    max_value = _nonlinear_binary
    min_value = _nonlinear_binary
    power = _nonlinear_binary

    zero(self, expr):
        raise Exception('Expressions containing zero are not allowed in IsMultiLinear')


def isMultiLinear(expr, arguments=None):
    if isinstance(expr, Form):
        check = IsMultiLinear(arguments if arguments is not None else expr.arguments())
        return all(r == 1 for r in map_expr_dags(check, (i.integrand() for i in expr.integrals())))
    elif isinstance(expr, Expr):
        check = IsMultiLinear(arguments if arguments is not None else extract_arguments(expr))
        return map_expr_dag(check, expr) == 1
    else:
        raise ArgumentError('Argument "expr" of isMultiLinear must be either of type Expr of of type Form.')
