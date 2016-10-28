from __future__ import division, print_function, unicode_literals

import copy, io, sys


# FileWriter
# ----------

class FileWriter:
    def __init__(self, fileName):
        self.file = open(fileName, "wt")

    def emit(self, src):
        print(src, file=self.file)

    def close(self):
        self.file.close()



# StringWriter
# ------------

class StringWriter:
    def __init__(self):
        self.file = io.StringIO()

    def emit(self, src):
        print(src, file=self.file)

    def getvalue(self):
        return self.file.getvalue()

    def close(self):
        self.file.close()



# ListWriter
# ----------

class ListWriter:
    def __init__(self):
        self.lines = []

    def emit(self, src):
        self.lines.append(src)

    def close(self):
        pass




# Block
# -----

class Block:
    def __init__(self):
        self.content = []

    def append(self, *objs):
        for obj in objs:
            if isinstance(obj, (list, set, tuple)):
                self.content += [o for o in obj]
            else:
                self.content.append(obj)



# Statement
# ---------

class Statement:
    def __init__(self):
        pass



# Expression
# ----------

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
            return Application(get(index), args=(self,))
        else:
            if isinstance(index, tuple):
                result = self
                for i in index:
                    result = Application(BracketOperator(), args=(result, makeExpression(i)))
                return result
            else:
                return Application(BracketOperator(), args=(self, makeExpression(index)))



# makeExpression
# --------------

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



# Application
# -----------

class Application(Expression):
    def __init__(self, function, args=None):
        Expression.__init__(self)
        self.function = function
        self.args = tuple(args)

    def __hash__(self):
        return hash((self.cppType, self.function, self.args))



# ConstantExpression
# ------------------

class ConstantExpression(Expression):
    def __init__(self, cppType, value):
        Expression.__init__(self, cppType)
        self.value = value

    def __hash__(self):
        return hash((self.cppType, self.value))



# ConstructExpression
# -------------------

class ConstructExpression(Expression):
    def __init__(self, cppType, args=None):
        Expression.__init__(self, cppType)
        self.args = None if args is None else [makeExpression(arg) for arg in args]

    def __hash__(self):
        return hash((self.cppType, self.args))



# InitializerList
# ---------------

class InitializerList(Expression):
    def __init__(self, *args):
        Expression.__init__(self)
        self.args = tuple(args)



# LambdaExpression
# ----------------

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



# Variable
# --------

class Variable(Expression):
  def __init__(self, cppType, name):
      Expression.__init__(self, cppType)
      self.name = name

  def __hash__(self):
      return hash((self.cppType, self.name))



# NameSpace
# ---------

class NameSpace(Block):
    def __init__(self, name=None):
        Block.__init__(self)
        self.name = name



# BuiltInFunction
# ---------------

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
        return Application(self, args=[makeExpression(arg) for arg in args])



# Operator
# --------

class Operator:
    def __init__(self, name, numArgs):
        self.name = name
        self.numArgs = numArgs



# PrefixUnaryOperator
# -------------------

class PrefixUnaryOperator(Operator):
    def __init__(self, name):
        Operator.__init__(self, name, 1)



# PostfixUnaryOperator
# -------------------

class PostfixUnaryOperator(Operator):
    def __init__(self, name):
        Operator.__init__(self, name, 1)



# BinaryOperator
# --------------

class BinaryOperator(Operator):
    def __init__(self, name):
        Operator.__init__(self, name, 2)



# BracketOperator
# ---------------

class BracketOperator(Operator):
    def __init__(self):
        Operator.__init__(self, '[]', 2)



# Function
# --------

class Function(Block):
    def __init__(self, cppType, name, targs=None, args=None, code=None):
        Block.__init__(self)
        self.cppType = cppType
        self.name = name
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() for a in args]
        if code is not None:
            self.append(code)



# Method
# ------

class Method(Block):
    def __init__(self, cppType, name, targs=None, args=None, code=None, static=False, const=False, volatile=False):
        Block.__init__(self)
        self.cppType = cppType
        self.name = name
        self.static = static
        self.resetQualifiers(const=const, volatile=volatile)
        if static and (const or volatile):
            raise Exception('Cannot cv-qualify static method.')
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() for a in args]
        if code is not None:
            self.append(code)

    def variant(self, typedName, const=False, volatile=False):
        method = copy.copy(self)
        method.typedName = typedName
        method.resetQualifiers(const=const, volatile=volatile)
        return method

    def resetQualifiers(self, const=False, volatile=False):
        self.qualifiers=[]
        if const:
            self.qualifiers.append('const')
        if volatile:
            self.qualifiers.append('volatile')



# Constructor
# -----------

class Constructor(Block):
    def __init__(self, targs=None, args=None, init=None, code=None):
        Block.__init__(self)
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() for a in args]
        self.init = None if init is None else [a.strip() for a in init]
        if code is not None:
            self.append(code)



# TypeAlias
# ---------

class TypeAlias:
    def __init__(self, name, typeName, targs=None):
        self.name = name
        self.typeName = typeName
        self.targs = None if targs is None else [a.strip() for a in targs]



# Declaration
# -----------

class Declaration:
    def __init__(self, obj, initializer=None, static=False, mutable=False):
        self.obj = obj
        self.initializer = initializer
        self.static = static
        self.mutable = mutable



class Using:
    def __init__(self, obj):
        self.obj = obj



# Class
# -----

class Class(Block):
    def __init__(self, name, targs=None, bases=None, final=False):
        Block.__init__(self)
        self.name = name
        self.targs = None
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.bases = None if bases is None else [a.strip() for a in bases]
        self.final = final
        self.type = 'class'

    def addMethod(self, typedName, targs=None, args=None, static=False, const=False, volatile=False, implemented=True):
        method = Method(typedName, targs=targs, args=args, static=static, const=const, volatile=volatile)
        if not implemented:
            method.append('DUNE_THROW(Dune::NotImplemented, "' + name + '::' + typedName + '");')
        append(method)
        return method if implemented else None



# Struct
# ------

class Struct(Class):
    def __init__(self, name, targs=None, bases=None, final=False):
        Class.__init__(self, name, targs=targs, bases=bases, final=final)
        self.type = 'struct'



# AccessModifier
# --------------

class AccessModifier:
    def __init__(self, modifier):
        if modifier not in ['private', 'protected', 'public']:
            raise Exception(modifier + ' is not a valid access modifier.')
        self.modifier = modifier



# EnumClass
# ---------

class EnumClass:
    def __init__(self, name, values, base=None):
        self.name = name
        self.base = None if base is None else base.strip()
        self.values = [v.strip() for v in values]



# ReturnStatement
# ---------------

class ReturnStatement(Statement):
    def __init__(self, expr=None):
        Statement.__init__(self)
        self.expression = None if expr is None else makeExpression(expr)



# short hand notation
# -------------------

def assign(left, right):
    return Application(BinaryOperator('='), args=(left, makeExpression(right)))


def construct(cppType, *args):
    return ConstructExpression(cppType, args if args else None)


def return_(expr=None):
    return ReturnStatement(expr)


def lambda_(args=None, capture=None, code=None):
    return LambdaExpression(args=args, capture=capture, code=code)


# built-in functions
# ------------------

atan = BuiltInFunction('cmath', 'X', 'atan', targs=['class X'], args=['conat X &x'])
atan2 = BuiltInFunction('cmath', 'X', 'atan2', targs=['class X'], args=['const X &x', 'const X &y'])
cos = BuiltInFunction('cmath', 'X', 'cos', targs=['class X'], args=['const X &x'])
pow_ = BuiltInFunction('cmath', 'X', 'pow', targs=['class X'], args=['const X &x', 'const X &y'])
sin = BuiltInFunction('cmath', 'X', 'sin', targs=['class X'], args=['const X &x'])
tan = BuiltInFunction('cmath', 'X', 'tan', targs=['class X'], args=['const X &x'])

def get(i):
    return BuiltInFunction('tuple', 'auto', 'get< ' + str(i) + ' >', targs=['class T'], args=['const T &arg'])

make_pair = BuiltInFunction('utility', 'std::pair< U, V >', 'make_pair', targs=['class U', 'class V'], args=['const U &left', 'const V &right'])



# SourceWriter
# ------------

class SourceWriter:
    def __init__(self, writer=None):
        if not writer:
            self.writer = StringWriter()
        elif self._isstring(writer):
            self.writer = FileWriter(writer)
        else:
            self.writer = writer
        self.blocks = []
        self.begin = True

    def close(self):
        if self.blocks:
            raise Exception("Open blocks left in source file.")
        self.writer.close()

    def emit(self, src, indent=0, context=None):
        if src is None:
            return
        elif isinstance(src, (list, set, tuple)):
            for s in src:
                self.emit(s, indent, context)
        elif isinstance(src, NameSpace):
            self.emit(None if self.begin else '')
            if src.name is not None:
                self.emit('namespace ' + src.name, indent)
            else:
                self.emit('namespace', indent)
            self.emit('{', indent)
            self.emit('', indent)
            self.emit(src.content, intent+1, src)
            self.emit('', indent)
            if src.name is not None:
                self.emit('} // namespace ' + src.name, indent)
            else:
                self.emit('} // anonymous namespace', indent)
        elif isinstance(src, Class):
            self.emit(None if self.begin else ['','',''], indent)
            self.emit(['// ' + src.name, '// ' + '-' * len(src.name), ''], indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            self.emit(src.type + ' ' + src.name + (' final' if src.final else ''))
            if src.bases:
                for i in range(0, len(src.bases)):
                    prefix = ': ' if i == 0 else '  '
                    postfix = ',' if i+1 < len(src.bases) else ''
                    self.emit(prefix + src.bases[i] + postfix, indent+1)
            self.emit('{', indent)
            self.emit(src.content, indent+1, src)
            self.emit('};', indent)
        elif isinstance(src, EnumClass):
            code = 'enum class ' + src.name
            if src.base is not None:
                code += ' : ' + src.base
            code += ' { ' + ', '.join(src.values) + ' };'
            self.emit(code, indent)
        elif isinstance(src, Function):
            self.emit(None if self.begin else '', indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            signature = self.typedName(src) + ' ('
            if src.args:
                signature += ' ' + ', '.join(src.args) + ' '
            signature += ')'
            self.emit(signature, indent)
            if src.content:
              self.emit('{', indent)
              self.emit(src.content, indent+1, src)
              self.emit('}', indent)
            else:
              self.emit('{}', indent)
        elif isinstance(src, Method):
            self.emit(None if self.begin else '', indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            signature = ('static ' if src.static else '') + self.typedName(src) + ' ('
            if src.args:
                signature += ' ' + ', '.join(src.args) + ' '
            signature += ')'
            if src.qualifiers:
                signature += ' ' + ' '.join(src.qualifiers)
            self.emit(signature, indent)
            if src.content:
              self.emit('{', indent)
              self.emit(src.content, indent+1, src)
              self.emit('}', indent)
            else:
              self.emit('{}', indent)
        elif isinstance(src, Constructor):
            if not isinstance(context, Class):
                raise Exception('Constructors can only occur in classes or structs')
            self.emit(None if self.begin else '', indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            signature = context.name + ' ('
            if src.args:
                signature += ' ' + ', '.join(src.args) + ' '
            signature += ')'
            self.emit(signature, indent)
            if src.init:
                for i in range(0,len(src.init)):
                    prefix = ': ' if i == 0 else '  '
                    postfix = ',' if i+1 < len(src.init) else ''
                    self.emit(prefix + src.init[i] + postfix, indent+1)
            if src.content:
              self.emit('{', indent)
              self.emit(src.content, indent+1, src)
              self.emit('}', indent)
            else:
              self.emit('{}', indent)
        elif isinstance(src, AccessModifier):
            if not isinstance(context, Class):
                raise Exception('Access modifiers can only occur in classes or structs')
            self.emit(None if self.begin else '')
            self.emit(src.modifier + ':', indent-1)
            self.begin = True
        elif isinstance(src, TypeAlias):
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
                self.emit('using ' + src.name + ' = ' + src.typeName + ';', indent)
            else:
                self.emit('typedef ' + src.typeName + ' ' + src.name + ';', indent)
        elif isinstance(src, Declaration):
            if not isinstance(src.obj, Variable):
                raise Exception('Only variables can be declared for now.')
            declaration = ('static ' if src.static else '') + ('mutable ' if src.mutable else '') + self.typedName(src.obj)
            if src.initializer is not None:
                expr = self.translateExpr(src.initializer)
                expr[len(expr)-1] += ';'
                self.emit(declaration + ' = ' + expr[0], indent)
                for e in expr[1:]:
                    if isinstance(e, tuple):
                        self.emit(e, indent+2, Function('auto', '<lambda>'))
                    else:
                        self.emit(e, indent+1)
            else:
                self.emit(declaration + ';', indent)
        elif isinstance(src, Using):
            if not isinstance(src.obj, BuiltInFunction):
                raise Exception('Only built-in functions can be used for now')
            if src.obj.namespace is None:
                function = src.obj.name
            else:
                function = src.obj.namespace + '::' + src.obj.name
            self.emit('using ' + function + ';')
        elif isinstance(src, Statement):
            if not isinstance(context, (Constructor, Function, Method)):
                raise Exception('Statements can only occur in constructors, functions and methods')
            if isinstance(src, Expression):
                expr = self.translateExpr(src)
                expr[len(expr)-1] += ';'
                self.emit(expr[0], indent)
                for e in expr[1:]:
                    if isinstance(e, tuple):
                        self.emit(e, indent+2, Function('auto', '<lambda>'))
                    else:
                        self.emit(e, indent+1)
            elif isinstance(src, ReturnStatement):
                if src.expression is not None:
                    expr = self.translateExpr(src.expression)
                    expr[len(expr)-1] += ';'
                    self.emit('return ' + expr[0], indent)
                    for e in expr[1:]:
                        if isinstance(e, tuple):
                            self.emit(e, indent+2, Function('auto', '<lambda>'))
                        else:
                            self.emit(e, indent+1)
                else:
                    self.emit('return;', indent)
            else:
                raise Exception('Unknown statement type.')
        elif self._isstring(src):
            src = src.rstrip()
            if src:
                self.writer.emit('  ' * (len(self.blocks) + indent) + src)
            else:
                self.writer.emit('')
            self.begin = False
        else:
            raise Exception("Unable to print " + repr(src) + ".")

    def translateExpr(self, expr):
        def join(lists, delimiter=''):
            left = lists[0]
            if len(lists) > 1:
                right = lists[1]
                n = len(left)
                return join([left[:n-1] + [left[n-1] + delimiter + right[0]] + right[1:]] + lists[2:], delimiter)
            else:
                return left

        if isinstance(expr, Application):
            args = [self.translateExpr(arg) for arg in expr.args]
            if isinstance(expr.function, Operator):
                if len(args) != expr.function.numArgs:
                    raise Exception('The operator' + expr.function.name + ' takes ' + expr.function.numArgs + ' arguments (' + str(len(args)) + ' given).')
                if isinstance(expr.function, PrefixUnaryOperator):
                    return join([[expr.function.name + '('], args[0], [')']])
                elif isinstance(expr.function, PostfixUnaryOperator):
                    return join([['('], args[0], [')' + expr.function.name]])
                elif isinstance(expr.function, BinaryOperator):
                    return join([['('], args[0], [') ' + expr.function.name + ' ('], args[1], [')']])
                elif isinstance(expr.function, BracketOperator):
                    return join([['('], args[0], [')[ '], args[1], [' ]']])
                else:
                    raise Exception('Unknown operator: ' + repr(expr.function))
            if isinstance(expr.function, BuiltInFunction):
                if expr.function.namespace is None:
                    function = expr.function.name
                else:
                    function = expr.function.namespace + '::' + expr.function.name
            elif isinstance(expr.function, (Function, Method)):
                function = expr.function.name
            else:
                function = expr.function
            if expr.args is not None:
                return join([[function + '( '], join(args, ', '), [' )']])
            else:
                return [function + '()']
        elif isinstance(expr, ConstantExpression):
            return [expr.value]
        elif isinstance(expr, ConstructExpression):
            if expr.args is not None:
                return join([[expr.cppType + '( '], join([self.translateExpr(arg) for arg in expr.args], ', '), [' )']])
            else:
                return [expr.cppType + '()']
        elif isinstance(expr, InitializerList):
            return join([['{ '], join([self.translateExpr(arg) for arg in expr.args], ', '), [' }']])
        elif isinstance(expr, LambdaExpression):
            capture = '' if expr.capture is None else ' ' + ', '.join(expr.capture) + ' '
            args = '' if expr.args is None else ' ' + ', '.join(expr.args) + ' '
            return ['[' + capture + '] (' + args + ') {', expr.code, '}']
        elif isinstance(expr, Variable):
            return [expr.name]
        else:
            return [expr.strip()]

    def typedName(self, obj):
        if obj.cppType.endswith('&') or obj.cppType.endswith('*'):
            return obj.cppType + obj.name
        else:
            return obj.cppType + ' ' + obj.name

    def pushBlock(self, block, name):
        self.blocks.append((block, name))
        self.begin = True

    def popBlock(self, block, name=None):
        (realBlock, realName) = self.blocks.pop()
        if realBlock != block or (name and name != realName):
            self.blocks.append((realBlock, realName))
            raise Exception("Trying to close " + realBlock + " " + realName + " as " + block + " " + name + ".");

    def getName(self,block):
        for b in reversed(self.blocks):
            if b[0] == block:
                return b[1]
        return None

    def openNameSpace(self, name):
        self.emit(None if self.begin else '')
        self.emit('namespace ' + name)
        self.emit('{')
        self.emit('')
        self.pushBlock('namespace', name)

    def closeNameSpace(self, name=None):
        self.emit(None if self.begin else '')
        self.popBlock('namespace', name)
        self.emit('} // namespace ' + name)

    def openClass(self, name, targs=None, bases=None):
        self.emit(None if self.begin else ['','',''])
        self.emit(['// ' + name, '// ' + '-' * len(name), ''])
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        self.emit('class ' + name)
        if bases:
            for i in range(0,len(bases)):
                prefix = '  : ' if i == 0 else '    '
                postfix = ',' if i+1 < len(bases) else ''
                self.emit(prefix + base + postfix)
        self.emit('{')
        self.pushBlock('class', name)

    def closeClass(self, name=None):
        self.popBlock('class', name)
        self.emit('};')

    def openStruct(self, name, targs=None, bases=None):
        self.emit(None if self.begin else ['','',''])
        self.emit(['// ' + name, '// ' + '-' * len(name), ''])
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        self.emit('struct ' + name)
        if bases:
            for i in range(0,len(bases)):
                prefix = '  : ' if i == 0 else '    '
                postfix = ',' if i+1 < len(bases) else ''
                self.emit(prefix + bases[i] + postfix)
        self.emit('{')
        self.pushBlock('struct', name)

    def closeStruct(self, name=None):
        self.popBlock('struct', name)
        self.emit('};')

    def section(self, section):
        self.emit(None if self.begin else '')
        (block, name) = self.blocks.pop()
        if block != 'class' and block != 'struct':
            self.blocks.push((block, name))
            raise Exception("Trying to declare " + section.strip() + " section outside class or struct.")
        self.emit(section.strip() + ':')
        self.blocks.append((block, name))
        self.begin = True

    def openConstMethod(self, typedName, targs=None, args=None, implemented=True):
        self.emit(None if self.begin else '')
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' ) const')
        else:
            self.emit(typedName + ' () const')
        self.emit('{')
        self.pushBlock('const method', typedName)
        if not implemented:
            struct = self.getName("struct")
            if struct:
                fullMethodName = struct + "::" + typedName
            else:
                fullMethodName = typedName
            self.emit('DUNE_THROW(Dune::NotImplemented, "'+fullMethodName+'");')

    def closeConstMethod(self, typedName=None):
        self.popBlock('const method', typedName)
        self.emit('}')

    def openMethod(self, typedName, targs=None, args=None):
        self.emit(None if self.begin else '')
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' )')
        else:
            self.emit(typedName + ' ()')
        self.emit('{')
        self.pushBlock('method', typedName)

    def closeMethod(self, typedName=None):
        self.popBlock('method', typedName)
        self.emit('}')

    def openFunction(self, typedName, targs=None, args=None):
        self.emit(None if self.begin else '')
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' )')
        else:
            self.emit(typedName + ' ()')
        self.emit('{')
        self.pushBlock('function', typedName)

    def closeFunction(self, typedName=None):
        self.popBlock('function', typedName)
        self.emit('}')

    def openRangeBasedFor(self, typedName, container):
        self.emit('for( ' + typedName.strip() + ' : ' + container.strip() + ' )')
        self.emit('{')
        self.pushBlock('range-based for', typedName)

    def closeRangeBasedFor(self, typedName=None):
        self.popBlock('range-based for', typedName)
        self.emit('}')

    def openPythonModule(self, moduleName):
        self.emit(None if self.begin else '')
        self.emit('PYBIND11_PLUGIN( ' + moduleName.strip() + ' )')
        self.emit('{')
        self.pushBlock('pybind11 module', moduleName)
        self.emit('pybind11::module module( "' + moduleName.strip() + '" );')

    def closePythonModule(self, moduleName=None):
        self.emit('')
        self.emit('return module.ptr();')
        self.popBlock('pybind11 module', moduleName)
        self.emit('}')

    def typedef(self, typeName, typeAlias, targs=None):
        self.emit(TypeAlias(typeAlias, typeName, targs=targs))

    def _isstring(self, obj):
        if sys.version_info.major == 2:
            return isinstance(obj, basestring)
        else:
            return isinstance(obj, str)

    @staticmethod
    def cpp_fields(field):
        if field=="complex":
            return "std::complex<double>"
        else:
            return "double"
