##############################
Information for C++ Developers
##############################

.. todo:: add something on compilerflag setting (see also `algorithm` section , e.g., explain

```
import dune.generator as generator
generator.addToFlags("-DWANT_CACHED_COMM_MANAGER=0",noChecks=True)
algorithm(...)
generator.setFlags("-g -Wfatal-errors",noChecks=True)
algorithm(...)
generator.reset()
```

.. todo:: mention use of `ccache` and `gdb`

.. todo:: mention `rmgenerated` script

.. todo:: mention fixes to `config.opts` e.g.

```
-DCMAKE_CXX_COMPILER=/usr/lib/ccache/g++
-DCMAKE_POSITION_INDEPENDENT_CODE=TRUE
-DALLOW_CXXFLAGS_OVERWRITE=ON
-DDUNE_PYTHON_INSTALL_EDITABLE=TRUE
-DADDITIONAL_PIP_PARAMS=\"-upgrade\"
```
