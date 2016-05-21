import sympy
from sympy.printing.ccode import CCodePrinter

# EXTENDING sympy.ccode
# make the sympy code generater use [i][j] notation instead of [i+n*j] as is used by defauls
class CCodePrinter(CCodePrinter):
    def _print_Indexed(self, expr):
        output = self._print(expr.base.label) + ''.join([ '[' + self._print(x) + ']' for x in expr.indices ])
        return output
def ccode(expr, assign_to=None, **settings):
    return CCodePrinter(settings).doprint(expr, assign_to)
