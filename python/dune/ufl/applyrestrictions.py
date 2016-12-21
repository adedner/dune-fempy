from ufl.corealg.multifunction import MultiFunction
from ufl.corealg.map_dag import map_expr_dag

class RestrictionPropagator(MultiFunction):
    def __init__(self, side=None):
        MultiFunction.__init__(self)
        self.side = side

    def _apply_restriction(self, expr):
        return expr if self.side is None else expr(self.side)

    def restricted(self, expr):
        # We allow for multiple restriction:
        # After restriction, the quantity is single-valued and cannot be restricted further.
        return map_expr_dag(RestrictionPropagator(expr.side()), expr.ufl_operands[0])

    operator = MultiFunction.reuse_if_untouched
    terminal = MultiFunction.reuse_if_untouched

    argument = _apply_restriction
    coefficient = _apply_restriction
    grad = _apply_restriction

    def variable(self, expr, op):
        return op

    def facet_normal(self, expr):
        return -expr('+') if self.side == '-' else expr('+')

    cell_volume = _apply_restriction


def applyRestrictions(expression):
    return map_expr_dag(RestrictionPropagator(), expression)