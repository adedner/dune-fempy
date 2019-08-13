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
renderView1.CenterOfRotation = [1.25, 1.25, 1.25]
renderView1.StereoType = 0
renderView1.CameraPosition = [4.234430847948241, -0.9444202633645626, -4.475037094059201]
renderView1.CameraFocalPoint = [0.5732842083224422, 1.7475819248138815, 2.548144673766888]
renderView1.CameraViewUp = [0.8982300684540759, 0.11454281719451433, 0.42433793980049694]
renderView1.CameraParallelScale = 2.165063509461097
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
spiral000 = XMLUnstructuredGridReader(FileName=['/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00011.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00012.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00013.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00014.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00015.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00016.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00017.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00018.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00019.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00020.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00021.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00022.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00023.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00024.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00025.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00026.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00027.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00028.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00029.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00030.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00031.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00032.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00033.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00034.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00035.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00036.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00037.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00038.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00039.vtu', '/home/dedner/DUNEPY3.core/dune-fempy/demos/spiral00040.vtu'])
spiral000.CellArrayStatus = ['estimate', 'levels', 'pDegree']
spiral000.PointArrayStatus = ['u', 'v']

# create a new 'Clip'
clip1 = Clip(Input=spiral000)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'u']
clip1.Value = 0.5050587803125381

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [1.25, 1.25, 1.25]
clip1.ClipType.Normal = [0.0, 1.0, 1.0]

# create a new 'Contour'
contour1 = Contour(Input=spiral000)
contour1.ContourBy = ['POINTS', 'u']
contour1.Isosurfaces = [0.5050587803125381]
contour1.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'levels'
levelsLUT = GetColorTransferFunction('levels')
levelsLUT.RGBPoints = [1.0, 0.231373, 0.298039, 0.752941, 2.5, 0.865003, 0.865003, 0.865003, 4.0, 0.705882, 0.0156863, 0.14902]
levelsLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'levels'
levelsPWF = GetOpacityTransferFunction('levels')
levelsPWF.Points = [1.0, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]
levelsPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour1
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.ColorArrayName = ['POINTS', '']
contour1Display.EdgeColor = [1.0, 1.0, 1.0]
contour1Display.GlyphType = 'Arrow'

# show data from clip1
clip1Display = Show(clip1, renderView1)
# trace defaults for the display properties.
clip1Display.Representation = 'Surface With Edges'
clip1Display.ColorArrayName = ['CELLS', 'levels']
clip1Display.LookupTable = levelsLUT
clip1Display.EdgeColor = [1.0, 1.0, 1.0]
clip1Display.GlyphType = 'Arrow'
clip1Display.ScalarOpacityUnitDistance = 0.12003354450516343
