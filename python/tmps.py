# trace generated using paraview version 5.10.0-RC2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
medium_bandsvtu = XMLUnstructuredGridReader(registrationName='medium_bands.vtu', FileName=['/home/remi/EmpiricalPseudopotential/build/medium_bands.vtu'])
medium_bandsvtu.PointArrayStatus = ['band_0', 'band_1', 'band_2', 'band_3', 'band_4', 'band_5', 'band_6', 'band_7', 'band_8', 'band_9', 'band_10', 'band_11']

# Properties modified on medium_bandsvtu
medium_bandsvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1_1 = CreateView('RenderView')
renderView1_1.AxesGrid = 'GridAxes3DActor'
renderView1_1.StereoType = 'Crystal Eyes'
renderView1_1.CameraFocalDisk = 1.0
renderView1_1.BackEnd = 'OSPRay raycaster'
renderView1_1.OSPRayMaterialLibrary = materialLibrary1

# show data in view
medium_bandsvtuDisplay = Show(medium_bandsvtu, renderView1_1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
medium_bandsvtuDisplay.Representation = 'Surface'
medium_bandsvtuDisplay.ColorArrayName = [None, '']
medium_bandsvtuDisplay.SelectTCoordArray = 'None'
medium_bandsvtuDisplay.SelectNormalArray = 'None'
medium_bandsvtuDisplay.SelectTangentArray = 'None'
medium_bandsvtuDisplay.OSPRayScaleArray = 'band_0'
medium_bandsvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
medium_bandsvtuDisplay.SelectOrientationVectors = 'None'
medium_bandsvtuDisplay.ScaleFactor = 0.2
medium_bandsvtuDisplay.SelectScaleArray = 'None'
medium_bandsvtuDisplay.GlyphType = 'Arrow'
medium_bandsvtuDisplay.GlyphTableIndexArray = 'None'
medium_bandsvtuDisplay.GaussianRadius = 0.01
medium_bandsvtuDisplay.SetScaleArray = ['POINTS', 'band_0']
medium_bandsvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
medium_bandsvtuDisplay.OpacityArray = ['POINTS', 'band_0']
medium_bandsvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
medium_bandsvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
medium_bandsvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
medium_bandsvtuDisplay.ScalarOpacityUnitDistance = 0.05082387603360301
medium_bandsvtuDisplay.OpacityArrayName = ['POINTS', 'band_0']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
medium_bandsvtuDisplay.ScaleTransferFunction.Points = [-12.58912229538956, 0.0, 0.5, 0.0, -8.17940466200778, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
medium_bandsvtuDisplay.OpacityTransferFunction.Points = [-12.58912229538956, 0.0, 0.5, 0.0, -8.17940466200778, 1.0, 0.5, 0.0]

# add view to a layout so it's visible in UI
AssignViewToLayout(view=renderView1_1, layout=None, hint=0)

# reset view to fit data
renderView1_1.ResetCamera(False)

# update the view to ensure updated data information
renderView1_1.Update()

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=medium_bandsvtu)
contour1.ContourBy = ['POINTS', 'band_0']
contour1.Isosurfaces = [-10.38426347869867]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1_1, 'GeometryRepresentation')

# get color transfer function/color map for 'band_0'
band_0LUT = GetColorTransferFunction('band_0')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'band_0']
contour1Display.LookupTable = band_0LUT
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = 'band_0'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.14208833774294388
contour1Display.SelectScaleArray = 'band_0'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'band_0'
contour1Display.GaussianRadius = 0.007104416887147194
contour1Display.SetScaleArray = ['POINTS', 'band_0']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'band_0']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-10.38426347869867, 0.0, 0.5, 0.0, -10.382309913635254, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-10.38426347869867, 0.0, 0.5, 0.0, -10.382309913635254, 1.0, 0.5, 0.0]

# hide data in view
Hide(medium_bandsvtu, renderView1_1)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1_1, True)

# update the view to ensure updated data information
renderView1_1.Update()

# get opacity transfer function/opacity map for 'band_0'
band_0PWF = GetOpacityTransferFunction('band_0')

# reset view to fit data
renderView1_1.ResetCamera(False)

# create new layout object 'Layout #2'
layout2 = CreateLayout(name='Layout #2')

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraFocalDisk = 1.0
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView2, layout=layout2, hint=0)

# set active view
SetActiveView(renderView1_1)

# destroy renderView1_1
Delete(renderView1_1)
del renderView1_1

# get layout
layout1 = GetLayoutByName("Layout #1")

RemoveLayout(layout1)

# set active source
SetActiveSource(contour1)

# show data in view
contour1Display = Show(contour1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'band_0']
contour1Display.LookupTable = band_0LUT
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = 'band_0'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.14208833774294388
contour1Display.SelectScaleArray = 'band_0'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'band_0'
contour1Display.GaussianRadius = 0.007104416887147194
contour1Display.SetScaleArray = ['POINTS', 'band_0']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'band_0']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [-10.38426347869867, 0.0, 0.5, 0.0, -10.382309913635254, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [-10.38426347869867, 0.0, 0.5, 0.0, -10.382309913635254, 1.0, 0.5, 0.0]

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView2, True)

# reset view to fit data
renderView2.ResetCamera(False)

# set scalar coloring
ColorBy(contour1Display, ('POINTS', 'band_1'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(band_0LUT, renderView2)

# rescale color and/or opacity maps used to include current data range
contour1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView2, True)

# get color transfer function/color map for 'band_1'
band_1LUT = GetColorTransferFunction('band_1')

# get opacity transfer function/opacity map for 'band_1'
band_1PWF = GetOpacityTransferFunction('band_1')

# rescale color and/or opacity maps used to exactly fit the current data range
contour1Display.RescaleTransferFunctionToDataRange(False, True)

# get animation track
contour1ContourValuesTrack = GetAnimationTrack('ContourValues', index=-1, proxy=contour1)

# create keyframes for this animation track

# create a key frame
keyFrame18356 = CompositeKeyFrame()
keyFrame18356.KeyValues = [-12.58912229538956]

# create a key frame
keyFrame18357 = CompositeKeyFrame()
keyFrame18357.KeyTime = 1.0
keyFrame18357.KeyValues = [-8.17940466200778]

# initialize the animation track
contour1ContourValuesTrack.KeyFrames = [keyFrame18356, keyFrame18357]

# get camera animation track for the view
cameraAnimationCue1 = GetCameraTrack(view=renderView2)

# create keyframes for this animation track

# create a key frame
keyFrame18363 = CameraKeyFrame()
keyFrame18363.Position = [2.1892280820889898e-05, -9.51454941421348e-05, 4.753830176434425]
keyFrame18363.FocalPoint = [2.1892280820889898e-05, -9.51454941421348e-05, -5.9556226088719466e-05]
keyFrame18363.ParallelScale = 1.2303972011298718
keyFrame18363.PositionPathPoints = [0.0, 0.0, 10.0, 7.818369630266645, 0.0, 6.234892711160037, 9.749363948615, 0.0, -2.2252608048500573, 4.338904848223074, 0.0, -9.009792394793298, -4.338821615061574, 0.0, -9.009811392219264, -9.749310421024916, 0.0, -2.225303491678726, -7.818353144918225, 0.0, 6.234858478981305]
keyFrame18363.FocalPathPoints = [2.18923e-05, -9.51455e-05, -5.95562e-05]
keyFrame18363.ClosedPositionPath = 1

# create a key frame
keyFrame18364 = CameraKeyFrame()
keyFrame18364.KeyTime = 1.0
keyFrame18364.Position = [2.1892280820889898e-05, -9.51454941421348e-05, 4.753830176434425]
keyFrame18364.FocalPoint = [2.1892280820889898e-05, -9.51454941421348e-05, -5.9556226088719466e-05]
keyFrame18364.ParallelScale = 1.2303972011298718

# initialize the animation track
cameraAnimationCue1.Mode = 'Path-based'
cameraAnimationCue1.KeyFrames = [keyFrame18363, keyFrame18364]

# hide color bar/color legend
contour1Display.SetScalarBarVisibility(renderView2, False)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
band_1LUT.ApplyPreset('Turbo', True)

# Hide orientation axes
renderView2.OrientationAxesVisibility = 0


SaveAnimation("tmp.avi", renderView1, ImageResolution=[1200, 1200],
    FontScaling='Do not scale fonts',
    FrameRate=5,
    Compression=True,
    Quality=2, 
    FrameWindow=[0, 10])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout2.SetSize(892, 592)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView2
renderView2.CameraPosition = [2.1892280820889898e-05, -9.51454941421348e-05, 4.753830176434425]
renderView2.CameraFocalPoint = [2.1892280820889898e-05, -9.51454941421348e-05, -5.9556226088719466e-05]
renderView2.CameraParallelScale = 1.2303972011298718

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).