.. role:: red

.. raw:: html

    <style> .red {color:red} </style>

.. raw:: html

  <h1> On the fly module construction </h1>

.. contents::
   :local:
   :depth: 2
   :backlinks: top


.. _constructionintro:

***********************************
Exporting C++ Classes 
***********************************

We start off by describing the concepts used to export a Dune C++ class to python.
Consider a class with some method we want to export 

.. code-block:: c++

  struct DuneClass 
  {
    DuneClass(...);
    const ReturnTypeOne& methodOne( ... ) const;
    ReturnTypeTwo methodTwo( ... );
    void methodThree( ... );
  };

..  admonition:: :red:`Definition`

  In the following we will to a C++ type as **regular** if they are
  automatically exported to python or imported from python using boost
  (e.g. *POD* and *std::string*) or by some other extension library
  (e.g. *Dune::FildVector*). 
  So *regular* types don't need any special treatment or wrapper classes to be used within python. 

First Step: Generating a wrapper class
=============================================

We start by specifying a class with the following structure:

.. code-block:: c++

  template <class Dune>
  struct WrapperClass
  {
    typedef Dune DuneType;
    WrapperClass(const std::shared_ptr<DuneType>& duneType) 
    : duneType_(duneType) 
    {}
    const std::shared_ptr<DuneType>& duneType() 
    { 
      return duneType_;
    }
    private:
    std::shared_ptr<DuneType> duneType_;
  };

We will use this to export *DuneClass* to python starting with

.. code-block:: c++

  typedef WrapperClass<DuneClass> Export;
  auto obj = boost::python::class_< Export, std::shared_ptr<Export> > ("PythonName", no_init);

.. _Section makeConstructor:

Second Step: adding constructors
===============================================

We use the option of boost python to define the python constructor out of
class:

.. code-block:: c++

  obj.def("__init__", make_constructor( makeWrapperclass ) );

with

.. code-block:: c++

  std::shared_ptr< WrapperClass<DuneClass> > makeWrapperClass(...)
  {
    return std::make_shared< WrapperClass<DuneClass> >
        ( std::make_shared<DuneClass>(...) );
  }

where the omitted parameter list corresponds to the parameter list of the
*DuneClass* constructor. We will consider a different list of parameters
below.

We will see later how this approach makes it easy to define a custom
constructor even in the code generation phase.

.. _Section ReturnValueConverter:

Third Step: return value conversion
=====================================

If the return values of a method to be exported is a **regular** type no
special treatment is necessary and exporting the method is straightforward
with boost python. Now assume that *RetrunValueOne* for example refers to
some other class with the Dune framework, e.g. an instance of a discrete function.
We make the following assumption:

..  admonition:: :red:`Assumption`

  All return values are either **regular** or *references (constant or non constant)*.
  That implies that any return value given as a copy can be handled by boost
  python directly or a direct export is provided in some extension module.

So taking the example from above *ReturnTypeOne* might refer to some class that
is not directly exported to python. We thus can not directly return the
corresponding reference stored in *DuneClass*. There are in fact two
reasons why this is problematic

#. The type is not known to boost python 
#. If the reference were to be exported to python and stored there, then
   this would result in a dangling reference once the actual instance of the
   wrapper class is not available anymore.

So we need to define a wrapper class for the return value that we can
export to python and we need to make sure that the reference to
*ReturnTypeOne* remains active long enough. We do this using wrapper classes of
the form

.. code-block:: c++

  template <class Container>
  struct ReturnOneWrapper
  {
    ReturnOneWrapper(const ReturnTypeOne &ret, const std::shared_ptr<Container>& container) 
    : ret_(ret), container_container) {}
    ReturnValue method(...);
    // further methods to export

    const ReturnTypeOne &get() const
    { return ret_; }
    ReturnTypeOne &get() const
    { return const_cast<ReturnTypeOne&>(ret_); }
    private:
    const ReturnTypeOne &ret_;
    shared_ptr<Container> container_;
  };

Depending on the actual class *ReturnTypeOne* the class will contain additional
methods to be exported to python. We will now need to do three things:

#. Depending on *ReturnTypeOne* use the correct wrapper class.
#. Convert the return value of *methodOne* 
#. Export the wrapper class with the corresponding methods to python (we
   always use *std::shared_ptr* as holder type.

To this end we use a template class that can be specialized on different return
values. The default implementation of this class is

.. code-block:: c++

  template <class R, class Container>
  struct ReturnValueConverter
  {
    typedef R type;
    static type convert(R& r, ...) 
    { return r; }
    static R &restore(type t)
    { return t; }
    static void registerType(const char* name) 
    {}
  };

This implementation of the class can be used for any **regular** type
returned as a copy but needs to be specialized over any type that needs a
wrapper class, e.g.,

.. code-block:: c++

  template <class Container>
  struct ReturnValueConverter< ReturnTypeOne, Container >
  {
    typedef std::shared_ptr< ReturnOneWrapper<Container> > type;
    static type convert(const ReturnTypeOne& r, const std::shared_ptr<Container> container) 
    { return type(r,container); }
    static ReturnTypeOne& get( type t )
    { return t->get(); }
    static void registerType(const char* name) 
    {
      auto obj = boost::python::cls_< type, std::shared_ptr<type> >(name, non_init);
      obj.def("method", &type::method);
    }
  };

Forth Step: converting parameters
===============================================

The *ReturnValueConverter* can also be used to convert back from the type
exported to python to the underlying C++ reference type (remember that we
only use the *ReturnValueConverter* for reference types). To this purpose
the structure contains a *get* method. In the example given above the
exported python class is *GridFunctionDiscrete*. This class stores a
reference the underlying Dune class and this can be accessed using the
*get* which makes the implementation on the *ReturnValueConverter*
straightforward. 

Fifth Step: adding methods
===============================================

If we assume that the wrapper and required *ReturnValueConverter* implementations
are available we can now export the *methodOne* on our *WrapperOne*:

.. code-block:: c++

  typedef ReturnValueConverter< ReturnTypeOne, DuneClass > Converter;
  obj.def("methodOne", &methodOneWrapper );
  Convert::registerType( "methodOne_ret" );

with (assuming here *DuneClass::methodOne* takes one parameter)

.. code-block:: c++

  typename Converter::typer methodOneWrapper( WrapperClass<DuneClass> *self, const Converter::type param )
  {
    return Converter::convert( self->duneType()->methodOne( Coverter::get(param) ), self->duneType() );
  }

In many cases that is all that is required and therefore we provide a
helper function to export. So to export *DuneClass::methodOne*, do the return value conversion, and
registration with boost python. So the above code would be replaced with the following line:

.. code-block:: c++

  defFromDune(obj, "methodOne", &DuneClass::methodOne );

.. note::

  The parameter list of the exported method will be identical to the parameter
  list of the method on the originsal class *DuneClass*. If the parameter list has to be
  different (e.g. some conversion is required) then the short hand version can not be
  used. To avoid having to use the *ReturnValueConverter* explicitely a hybrid approach
  to the above could be useful:

  .. code-block:: c++

    const ReturnTypeOne &methodOneWrapper( DuneClass *self, MyParameterList... )
    {
      // convert MyParameterList to DuneParameterList 
      return self->methodOne(DuneParameterList);
    }

    defFromDune(obj, "methodOne", &methodOneWrapper );

  Note that here *self* is a pointer to the *DuneClass* type and the return
  type is the original type from *methodOne*.

.. _Summary:

Summary
=======================

In the above we distinguished three different C++ types

#. Types that boost python provides an export
   mechanism by default or other classes that provide boost python exports
   that can be used directly, i.e., without wrapper classes
   (we referred to these as **regular**).
#. Internal Dune classes that need to be wrapped before they can be
   exported to python. These classes can not be initialized directly from
   within python, i.e., they have no constructor. Thus instances of these
   classes will only appear as return values of functions or methods.
   The wrapper class internally stores the reference to an instance of the actual 
   Dune class and a shared pointer to a *container* which in most cases
   will be the variable that constructed this reference. This provides
   reference counting to avoid dangling references. The conversion from the
   Dune class to the exported wrapper is done using the
   *ReturnValueConverter* structure as described above. The wrapper should be derived
   from:

  .. code-block:: c++

    template <class Object, class Container>
    struct ReferenceWrapperClass
    {
      ReferenceWrapperClass(const Object &ret, const std::shared_ptr<Container>& container);
      const Object &get() const;
      Object &get();
    };

#. Finally we need to distinguish classes that need to be constructed from
   within python. These use wrappers for the Dune structures which are
   constructed from a shared pointer of the Dune class. The following wrapper should be
   used:

   .. code-block:: c++

    template <class Dune>
    struct WrapperClass
    {
      typedef Dune DuneType;
      WrapperClass(const std::shared_ptr<DuneType>& duneType);
      std::shared_ptr<DuneType> get() const;
    };

  This third category of types will be the focus of the following discussion. 

****************************************************************
Exporting Different Implementations of an Interface Class
****************************************************************

In the following we describe the approach taken for auto generating the
wrapper code for classes realizing a common interface. These will be types belonging to the third 
category described in the Summary_ above, i.e.,, Dune grid classes. We will
use the grid class as a show case in the following.
At first We will assume that we have a wrapper class that provides the common interface, e.g.,

.. code-block:: c++

  template <class Impl>
  struct GridWrapper
  {
    GridWrapper(const std::shared_ptr<InterfaceWrapper<Impl>> &);
    const ReturnTypeOne& methodOne( ... ) const;
    ReturnTypeTwo methodTwo( ... );
    void methodThree( ... );
  };

that is to be exported to python. This is implemented and exported using
boost python as described above, i.e., using the *ReturnValueConverter* and
the *make_constructor* mechanism of boost python. In the following we first
describe how the user can add a new implementation of this interface and
how just-in-time compilation of the wrapper is achieved. After that we will
explain how extensions of that interface can be provided by the user.

Database approach
===================

To provide the information required to export an implementation of the
wrapper class we use python dictonary containing information to generate the
typedef for the implementation and the required include statements needed
to compile the wrapper. For a grid examples given above files containing the required dictonaries 
would be located in the directory *python/database/grid*. An example would
be

.. code-block:: python

  "ALUCubeGrid" :
  {
    "type"      : "Dune::ALUGrid< $(dimgrid), $(dimworld), Dune::cube, Dune::nonconforming >",
    "default"   : [ "dimworld=$(dimgrid)" ],
    "checks"    : [ "$(dimgrid)==$(dimworld) or $(dimgrid)+1==$(dimworld)",
                    "2<=$(dimworld) and $(dimworld)<=3",
                    "2<=$(dimgrid) and $(dimgrid)<=3" ],
    "include"   : [ "dune/alugrid/grid.hh", "dune/alugrid/dgf.hh" ]
  }

The first line provides the information on how to construct the correct
typedef for the Dune class:

.. code-block:: c++

  template< int dimgrid, int dimworld,  ALU3dGridElementType elType>
  class ALU3dGrid<dimgrid,dimworld,Dune::cube,Dune::nonconforming>;

The last line contains information about which include statements are
required to compile the wrapper, i.e., 

.. code-block:: c++

   #include <dune/alugrid/grid.hh> 
   #include <dune/alugrid/dgf.hh>

The remaining two lines provide default values and tests for free
parameters *dimgrid* and *dimworld*.

On the fly code generation is done with the python module *dune.generator*. 

.. autoclass:: dune.generator.generator.Generator
   :members:
   :special-members: __init__

In the above example the extnesion module for our grid implementation would
be generated by calling

.. code-block:: python

     aluModule = Generator("grid").getModule("ALUCubeGrid", dimgrid=2, dimworld=3)

The following header file is generated in this example: 

.. code-block:: c++

  #include <dune/alugrid/grid.hh> 
  #include <dune/alugrid/dgf.hh>

  #include <dune/fempy/pygrid.hh>

  typedef ALU3dGrid<2,3,Dune::cube,Dune::nonconforming> DuneType;

  BOOST_PYTHON_MODULE( grid_a1fb60a2a2b9998d3d2288f8c7c39128 ) { registerDuneGrid< DuneType >(); }

The first two lines are the two includes given in the dictonary and the
typedef is constructed using the dictonary together with the named
parameters in the call to *getModule*. The third include is for the file
containing the wrapper class and the actual boost python exports statements in the free function 
*registerDuneGrid*. Note that the grid identifier used here is taken from
the constructor argument used to set up the *Generator* class. The
*registerDuneGrid* function can refer to the correct implementation using
the *DuneType* typedef, so the following class will be exported;

.. code-block:: c++
  
  auto obj = boost.python.class_< GridWrapper<DuneType>, std::shared_ptr<GridWrapper<DuneType> >(name,no_init);

The actual library for the extension mdoule will be generated in *python.dune.generated* within the
cmake build directory tree. The name of the library is taken by hashing the
Dune type of the implementation class. 

Summary
------------

The 
:py:class:`.generator.Generator`
class is used to perform on the fly
generation and import of extension modules for one given realization of a C++ interface class *InterfaceClass*. 
The details of each realization is provided through dictonaries contained
in files within the directory *python/database/interfaceName* where
*interfaceName* is the identifier used for the *InterfaceClass*. The
dictionary contains the means for constructing the correct type (*DuneType*
in the following) and a list of required include files. Calling the method
:py:meth:`.Generator.getModule` 
with the identifier used in
the dictionary for the desired implementation and the name parameters
required to fix *DuneType* results in a single header 
*python/dune/generated/generated_module.hh*. This file includes all headers
provided in the dictionary plus a header *pyinterfaceName.hh* containing the
definition of the wrapper code and the boost python export. The extension library is build 
using the *cmake* build system with the *generate_module* target. The resulting library
is renamed and imported into the python environment. 

.. note:: 
  The extension module will only be build if it does not exist
  already, i.e., no additional checks are performed to determine if
  the dependencies for this module have changed. 

Extending on the fly generated modules
==========================================

The method :py:meth:`.Generator.getModule` returns the python module after
it has been compiled and imported. At this point it is straightforward to
add attributes directly within python. This approach is used for example
for the grid modules (see for example the implementation of the function 
:py:func:`.fem.grid.get`). Extending the wrapped class contained in the
extension module (e.g. :py:class:`.fem.grid.LeafGrid`) using python is just
as straightforward.

Extending the wrapped class on the C++ side, i.e., before the generation of
the extension module requires a bit more work:

First we note that by using the approach described in the 
:ref:`Section on return value conversion <Section ReturnValueConverter>`
the return values of the interface method does not have to be part of the
interface but can be changed by the implementation. For any method for
which the export makes use of the 
*defFromDune* method also described in the same section the parameters
passed to the method can also be changed directly within the implementation
class. By either deciding to use or not to use this functionality the
writer of the wrapper class and code for the boost python export can makes
a decision on whether the exact method parameters and return values are to
be part of the interface or not. Depending on this choice either or both
the parameters and return values can be changed by the interface
realization. 

Next we describe how new methods can be added to the wrapper class:

To this end a free standing function *bool addToPython(...)* is provided.
To add new methods for an implementation the developer has to provide an overload of 
this function. Using again the example of grid classes used previously:

.. code-block:: c++

  template< class... Args >
  bool addToPython (boost::python::class_< GridWrapper<MyGrid>, Args... > &cls )
  { cls.def(...); }

This overload should be implemented in one of the include files provided in
the dictionary entry of the implementation.
Note that the *ReturnValueConverter* and the *defFromDune* function can be
used to simplify the export. 

So adding new method methods is also straightforward.

Finally, it remains to discuss changing the constructor of the implementation class. Normally a 
standard constructor will be exported using *boost.python.make_constructor* as described in Section on 
:ref:`adding constructors <Section makeConstructor>`.
This will be done after the call to *exportToPython". So no constructer
will have been added before the call to this function and the developer of
the interface implementation can use *make_constructor* to add a different custom
constructor. To avoid having the default constructor also being added after
returning from the *addToPython* call the function should return false. If
the default constructor is to be used the function should return true.


*********************************
Example: Interfaces Provided 
*********************************

Grid Construction
=======================

Module
-------

.. automodule:: dune.fem.grid
   :members: get, leafGrid
.. automodule:: griddocu
   :noindex:
   :members:
   :undoc-members:

Grid class
-----------

.. autoclass:: griddocu.LeafGrid
   :members:
   :undoc-members:
   :inherited-members:
   :special-members: __init__

VTK Output 
-------------

.. autoclass:: griddocu.VTKOutput
   :members:
   :undoc-members:
   :inherited-members:
   :special-members: __init__

Scheme Construction
=========================

Module
---------

.. automodule:: dune.fem.scheme
   :members: get, scheme
.. automodule:: schemedocu
   :noindex:
   :members:
   :undoc-members:

Scheme class
--------------

.. autoclass:: schemedocu.Scheme
   :members:
   :undoc-members:
   :inherited-members:
   :special-members: __init__

`Scheme wrapper class`_
----------------------------------------------
.. _Scheme wrapper class: file:../../doxygen/html/struct_scheme.html
