from dune.source.builtin import get, make_shared
from dune.source.cplusplus import AccessModifier, Constructor, Declaration, Function, Method, NameSpace, Struct, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import UnformattedExpression
from dune.source.cplusplus import SourceWriter
from dune.source.cplusplus import assign, construct, dereference, lambda_, nullptr, return_, this
from dune.models.elliptic.model import EllipticModel
from dune.ufl.codegen import CodeGenerator
from dune.source.fem import declareFunctionSpace

class SplitDomainModel(EllipticModel):
    def __init__(self, dimRange, signature):
        EllipticModel.__init__(self, dimRange, signature)
        self.alpha = ''
        self.linAlpha = ''
        self.initIntersection += '\n      initIntersectionImpl( intersection, std::make_index_sequence< numCoefficients >() );'
        self.extraMethods = self.intersectionImpl()
        #self.initCoefficients = 'std::ignore = std::make_tuple( (gridCheck(std::get< i >( coefficients_ ), i), 1) ... );'
        self.extraMethods.append(self.initCoefficients())
        self.extraMethods.append(self.gridCheck())
        self.vars += 'mutable std::vector< int > coeffInitialized_ = std::vector< int >( numCoefficients );'

    def coefficient(self, idx, x):
        coefficient = []
        for t, n in (('RangeType', 'evaluate'), ('JacobianRangeType', 'jacobian'), ('HessianRangeType', 'hessian')):
            result = Variable('typename std::tuple_element_t< ' + str(idx) + ', CoefficientFunctionSpaceTupleType >::' + t, 'result')
            code = [Declaration(result)]
            code.append('if (coeffInitialized_[' + str(idx) + '] == 1)')
            code.append('  std::get< ' + str(idx) + ' >( coefficients_ ).' + n + '( x, ' + result.name + ' );')
            code.append('else if (coeffInitialized_[' + str(idx) + '] == -1)')
            code.append('{')
            code.append('  using Dune::Fem::coordinate;')
            code.append('  auto xc = coordinate(x);')
            code.append('  auto xinter = intersection_->geometryInInside().local(xc);')
            code.append('  auto xnb = intersection_->geometryInOutside().global(xinter);')
            code.append('  std::get< ' + str(idx) + ' >( coefficients_ ).' + n + '( xnb, ' + result.name + ' );')
            code.append('}')
            code.append('else')
            code.append('{')
            code.append('  std::cout << "coefficient not initialized!" << std::endl;')
            code.append('  abort();')
            code.append('}')
            code.append('return result;')
            coefficient += [lambda_(capture=[this], args=['auto x'], code=code)(x)]
        return coefficient

    def intersectionImpl(self):
        output = []
        code1 = ['if (intersection.neighbor())']
        code1.append('{')
        code1.append('  const EntityType &nb = intersection.outside();')
        code1.append('  std::ignore = std::make_tuple( (initCoeffOnNb(nb, std::get< i >( coefficients_ ), i), 1) ... );')
        code1.append('}')
        output.append(Method('void', 'initIntersectionImpl', targs=['std::size_t... i'],
                        args=['const IntersectionType &intersection', 'std::index_sequence< i... >'], code=code1, const=True))
        code2 = ['if (coeffInitialized_[i] != 1)']
        code2.append('{')
        code2.append('  coeff.init( nb );')
        code2.append('  coeffInitialized_[i] = -1;')
        code2.append('}')
        output.append(Method('void', 'initCoeffOnNb', targs=['class CoeffType'],
                      args=['const EntityType &nb', 'CoeffType &coeff', 'std::size_t i'], code=code2, const=True))
        return output

    def initCoefficients(self):
        code = ['std::ignore = std::make_tuple( (gridCheck(std::get< i >( coefficients_ ), i), 1) ... );']
        return Method('void', 'initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'], code=code, const=True)

    def gridCheck(self):
        code = ['const auto &gridPart = coeff.gridPart();']
        code.append('if (gridPart.contains( entity() ))')
        code.append('{')
        code.append('  coeffInitialized_[i] = 1;')
        code.append('  coeff.init( entity() );')
        code.append('}')
        code.append('else')
        code.append('  coeffInitialized_[i] = 0;')
        return Method('void', 'gridCheck', targs=['class CoeffType'], args=['CoeffType &coeff', 'std::size_t i'], code=code, const=True)

    def code(self, name='Model', targs=[]):
        constants_ = Variable('std::tuple< ' + ', '.join('std::shared_ptr< ' + c  + ' >' for c in self._constants) + ' >', 'constants_')
        coefficients_ = Variable('std::tuple< ' + ', '.join('Coefficient' + str(i) for i, c in enumerate(self._coefficients)) + ' >', 'coefficients_')
        entity_ = Variable('const EntityType *', 'entity_')

        code = Struct(name, targs=(['class GridPart'] + ['class Coefficient' + str(i) for i, c in enumerate(self._coefficients)] + targs))

        code.append(TypeAlias("GridPartType", "GridPart"))
        code.append(TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"))
        code.append(TypeAlias("IntersectionType", "typename GridPart::IntersectionType"))

        code.append(declareFunctionSpace("typename GridPartType::ctype", SourceWriter.cpp_fields(self.field), UnformattedExpression("int", "GridPartType::dimensionworld"), self.dimRange))
        code.append(Declaration(Variable("const int", "dimLocal"), initializer=UnformattedExpression("int", "GridPartType::dimension"), static=True))

        if self.hasConstants:
            code.append(TypeAlias("ConstantType", "typename std::tuple_element_t< i, " + constants_.cppType + " >::element_type", targs=["std::size_t i"]))
            code.append(Declaration(Variable("const std::size_t", "numConstants"), initializer=len(self._constants), static=True))

        if self.hasCoefficients:
            coefficientSpaces = ["Dune::Fem::FunctionSpace< DomainFieldType, " + SourceWriter.cpp_fields(c['field']) + ", dimDomain, " + str(c['dimRange']) + " >" for c in self._coefficients]
            code.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " + ", ".join(coefficientSpaces) + " >"))
            code.append('static const std::size_t numCoefficients = std::tuple_size< CoefficientFunctionSpaceTupleType >::value;')
            code.append(TypeAlias('CoefficientType', 'std::tuple_element_t< i, ' + coefficients_.cppType + ' >', targs=['std::size_t i']))

        if self.hasCoefficients:
            args = [Variable("const Coefficient" + str(i) + " &", "coefficient" + str(i)) for i, c in enumerate(self._coefficients)]
            init = ["coefficients_( " + ", ".join("coefficient" + str(i) for i, c in enumerate(self._coefficients)) + " )"]
            constructor = Constructor(args=args, init=init)
        else:
            constructor = Constructor()
        if self.hasConstants:
            constructor.append([assign(get(str(i))(constants_), make_shared(c)()) for i, c in enumerate(self._constants)])
        code.append(constructor)

        init = ['entity_ = &entity;']
        init += ['initCoefficients( std::make_index_sequence< numCoefficients >() );']
        init = [UnformattedBlock(init)] + self.init + [return_(True)]
        code.append(Method('bool', 'init', args=['const EntityType &entity'], code=init, const=True))

        code.append(Method('const EntityType &', 'entity', code=return_(dereference(entity_)), const=True))
        code.append(Method('std::string', 'name', const=True, code=return_(UnformattedExpression('const char *', '"' + name + '"'))))

        code.append(TypeAlias("BoundaryIdProviderType", "Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >"))
        code.append(Declaration(Variable("const bool", "symmetric"), initializer=self.symmetric, static=True))

        code.append(Method('void', 'source', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.source, const=True))
        code.append(Method('void', 'linSource', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.linSource, const=True))

        code.append(Method('void', 'diffusiveFlux', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.diffusiveFlux, const=True))
        code.append(Method('void', 'linDiffusiveFlux', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.linDiffusiveFlux, const=True))

        code.append(Method('void', 'fluxDivergence', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_d2u, self.arg_r], code=self.fluxDivergence, const=True))

        code.append(Method('void', 'alpha', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_r], code=self.alpha, const=True))
        code.append(Method('void', 'linAlpha', targs=['class Point'], args=[self.arg_ubar, self.arg_x, self.arg_u, self.arg_r], code=self.linAlpha, const=True))

        code.append(Method('bool', 'hasDirichletBoundary', const=True, code=return_(self.hasDirichletBoundary)))
        code.append(Method('bool', 'hasNeumanBoundary', const=True, code=return_(self.hasNeumanBoundary)))

        code.append(Method('bool', 'isDirichletIntersection', args=[self.arg_i, 'Dune::FieldVector< int, dimRange > &dirichletComponent'], code=self.isDirichletIntersection, const=True))
        code.append(Method('void', 'dirichlet', targs=['class Point'], args=[self.arg_bndId, self.arg_x, self.arg_r], code=self.dirichlet, const=True))

        code.append(Method('void', 'initIntersection', args=['const IntersectionType &intersection'], code=self.initIntersection, const=True))
        code.append(self.extraMethods)

        if self.hasConstants:
            code.append(Method("const ConstantType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_))), const=True))
            code.append(Method("ConstantType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_)))))

        if self.hasCoefficients:
            code.append(Method("const CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_)), const=True))
            code.append(Method("CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_))))

        code.append(AccessModifier("private"))
        code.append(Declaration(entity_, nullptr, mutable=True))
        code.append(Declaration(Variable('const IntersectionType *', 'intersection_'), nullptr, mutable=True))
        if self.hasConstants:
            code.append(Declaration(constants_, mutable=True))
        if self.hasCoefficients:
            code.append(Declaration(coefficients_, mutable=True))
        code.append(self.vars)
        return code

class SplitDomainCodeGenerator(CodeGenerator):
    def __init__(self, predefined, coefficients, tempVars):
        CodeGenerator.__init__(self, predefined, coefficients, tempVars)

    def coefficient(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            pass

        print('Warning: ' + ('Constant ' if expr.is_cellwise_constant() else 'Coefficient ') + str(expr) + ' not predefined.')
        idx = str(self._getNumber(expr))
        if expr.is_cellwise_constant():
            var = Variable('const ConstantsRangeType< ' + idx + ' >', 'cc' + idx)
            self.code.append(Declaration(var, 'constant< ' + idx + ' >()'))
        else:
            var = Variable('CoefficientRangeType< ' + idx + ' >', 'c' + idx)
            self.code.append(Declaration(var))
            self.code.append('if (coeffInitialized_[' + idx + '] == 1)')
            self.code.append('coefficient< ' + idx + ' >().evaluate( x, c' + idx + ' );')
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
            self.code.append('  std::cout << "coefficient not initialized!" << std::endl;')
            self.code.append('  abort();')
            self.code.append('}')
        return var
