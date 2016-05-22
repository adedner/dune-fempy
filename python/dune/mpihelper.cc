#include <config.h>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/fempy/pybind11/pybind11.h>

PYBIND11_PLUGIN( mpihelper )
{
  pybind11::module module( "mpihelper" );

  try
  {
    int argc = 0;
    char **argv = nullptr;
    Dune::MPIHelper::instance( argc, argv );

    // managed to initialize mpi
    typedef Dune::CollectiveCommunication< Dune::MPIHelper::MPICommunicator > Comm;

    pybind11::class_< Comm > cc( module, "CollectiveCommunication" );
    cc.def_property_readonly( "rank", &Comm::rank );
    cc.def_property_readonly( "size", &Comm::size );
    cc.def( "barrier", &Comm::barrier );

    module.attr( "comm" ) = pybind11::cast( Dune::MPIHelper::getCollectiveCommunication() );
  }
  catch ( const std::exception &e )
  {
    std::cout << e.what() << std::endl;
  }

  return module.ptr();
}
