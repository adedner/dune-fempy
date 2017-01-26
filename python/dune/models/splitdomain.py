from dune.source import Method
from dune.models.elliptic import EllipticModel
from dune.ufl.codegen import CodeGenerator

class SplitDomainModel(EllipticModel):
    def __init__(self, dimRange, signature):
        EllipticModel.__init__(self, dimRange, signature)
        self.alpha = ''
        self.linAlpha = ''
        self.initIntersection += '\n      initIntersectionImpl( intersection, std::make_index_sequence< numCoefficients >() );'
        self.extraMethods = self.intersectionImpl()
        self.initCoefficients = 'std::ignore = std::make_tuple( (gridCheck(std::get< i >( coefficients_ ), i), 1) ... );'
        self.privateMethods = self.gridCheck()
        self.vars += 'mutable std::vector< int > coeffInitialized_ = std::vector< int >( numCoefficients );'

    def intersectionImpl(self):
        output = []
        code1 = ['if (intersection.neighbor())']
        code1.append('{')
        code1.append('  const EntityType &nb = intersection.outside();')
        code1.append('  std::ignore = std::make_tuple( (initCoeffOnNb(nb, std::get< i >( coefficients_ ), i), 1) ... );')
        code1.append('}')
        output.append(Method('void initIntersectionImpl', targs=['std::size_t... i'],
                        args=['const IntersectionType &intersection', 'std::index_sequence< i... >'], code=code1, const=True))
        code2 = ['if (coeffInitialized_[i] != 1)']
        code2.append('{')
        code2.append('  coeff.init( nb );')
        code2.append('  coeffInitialized_[i] = -1;')
        code2.append('}')
        output.append(Method('void initCoeffOnNb', targs=['class CoeffType'],
                      args=['const EntityType &nb', 'CoeffType &coeff', 'std::size_t i'], code=code2, const=True))
        return output

    def gridCheck(self):
        code = ['const auto &gridPart = coeff.gridPart();']
        code.append('if (gridPart.contains( entity() ))')
        code.append('{')
        code.append('  coeffInitialized_[i] = 1;')
        code.append('  coeff.init( entity() );')
        code.append('}')
        code.append('else')
        code.append('  coeffInitialized_[i] = 0;')
        return Method('void gridCheck', targs=['class CoeffType'], args=['CoeffType &coeff', 'std::size_t i'], code=code, const=True)

class SplitDomainCodeGenerator(CodeGenerator):
    def __init__(self, predefined, coefficients, tempVars):
        CodeGenerator.__init__(self, predefined, coefficients, tempVars)

    def coefficient(self, expr):
        try:
            return self.predefined[expr]
        except KeyError:
            pass

        idx = str(self._getNumber(expr))
        if expr.is_cellwise_constant():
            init = 'ConstantsRangeType< ' + idx + ' > cc' + idx + ' = constant< ' + idx + ' >();'
            if not init in self.code:
                self.code.append(init)
            return 'cc' + idx
        else:
            init = 'CoefficientRangeType< ' + idx + ' > c' + idx + ';'
            if not init in self.code:
                self.code.append(init)
                self.code.append('if (coeffInitialized_[' + idx + '] == 1)')
                self.code.append('  coefficient< ' + idx + ' >().evaluate( x, c' + idx + ' );')
                self.code.append('else if (coeffInitialized_[' + idx + '] == -1)')
                self.code.append('{')
                self.code.append('  using Dune::Fem::coordinate;')
                self.code.append('  auto xc = coordinate(x);')
                self.code.append('  auto xinter = intersection_->geometryInInside().local(xc);')
                self.code.append('  auto xnb = intersection_->geometryInOutside().global(xinter);')
                self.code.append('  coefficient< ' + idx + ' >().evaluate( xnb, c' + idx + ' );')
                self.code.append('}')
                self.code.append('else')
                self.code.append('{')
                self.code.append('  std::cout << "coefficient in alpha method not initialized!" << std::endl;')
                self.code.append('  abort();')
                self.code.append('}')
            return 'c' + idx
