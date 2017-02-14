# coding: utf-8

# # Parameters [Notebook][2]
# Parameters can either be read from parameter files or directly added during the python session. A parameter can only be added once and not changed afterwards. The parameters can be read by any part of the C++ code (or python code for that matter. For indvidual classes the default parameters used can be replaced by passing in an additional dictonary during construction. Parameters defined in that dictonary will replace the global parameters which will be used as default value for any parameter missing in the user provided dictonary. In this example we use a simple [parameter file][1].
#
# [1]: _downloads/parameter
# [2]: _downloads/parameters.ipymn

# In[ ]:

from dune.fem import parameter
parameter.append( "parameter" )
parameter.append( {"hi": 12, "how are": 20, "you": "getting the hang of it?" } )
parameter.append( "hi", 11. )
parameter["okay"] = 12
parameter.append( "I", "am fine" )
print(parameter["really"])


# All parameters can be easily printed. They are sorted first by the local from which they were read (files or program code). Parameters proceeded by a `#` have not yet been read.

# In[ ]:

print(parameter)
