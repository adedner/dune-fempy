# state file generated using paraview version 5.0.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1375, 662]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [-0.011274993419647217, 0.13851499557495117, 0.0009349584579467773]
renderView1.StereoType = 0
renderView1.CameraPosition = [5.943631346939618, -2.6750183863733055, -0.9818136264017825]
renderView1.CameraFocalPoint = [-2.225821723798293, 1.1848288745315703, 0.3664054698302595]
renderView1.CameraViewUp = [0.44643227403377556, 0.8653981632530174, 0.22756151638476937]
renderView1.CameraParallelScale = 2.364423238022929
renderView1.Background = [1.0, 1.0, 1.0]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Unstructured Grid Reader'
mcfSoap000 = XMLUnstructuredGridReader(FileName=['/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00000.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00001.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00002.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00003.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00004.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00005.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00006.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00007.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00008.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00009.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00010.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00011.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00012.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00013.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00014.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00015.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00016.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00017.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00018.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00019.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00020.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00021.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00022.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00023.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00024.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00025.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00026.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00027.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00028.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00029.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00030.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00031.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00032.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00033.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00034.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00035.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00036.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00037.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00038.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00039.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00040.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00041.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00042.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00043.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00044.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00045.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00046.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00047.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00048.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00049.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00050.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00051.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00052.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00053.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00054.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00055.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00056.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00057.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00058.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00059.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00060.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00061.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00062.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00063.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00064.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00065.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00066.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00067.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00068.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00069.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00070.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00071.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00072.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00073.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00074.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00075.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00076.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00077.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00078.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00079.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00080.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00081.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00082.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00083.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/mcfSoap00084.vtu'])
mcfSoap000.PointArrayStatus = ['uh']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from mcfSoap000
mcfSoap000Display = Show(mcfSoap000, renderView1)
# trace defaults for the display properties.
mcfSoap000Display.ColorArrayName = ['POINTS', '']
mcfSoap000Display.EdgeColor = [1.0, 1.0, 1.0]
mcfSoap000Display.GlyphType = 'Arrow'
mcfSoap000Display.ScalarOpacityUnitDistance = 0.12765646545333648
