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

... todo:: mention use of `ccache` and `gdb`

... todo:: mention `rmgenerated` script
