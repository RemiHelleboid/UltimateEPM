# trace generated using paraview version 5.10.0-RC2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

# import the simple module from the paraview
from paraview.simple import *
# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import argparse


    # parse command line argument
parser = argparse.ArgumentParser(prog='PROG',description='Test')
parser.add_argument("-f",dest="file",required=True, type=str)
parser.add_argument("-b",dest="band",required=False, type=int, default=0)
parser.add_argument("-e",dest="iso_energy",required=False, type=float, default=0)
parser.add_argument("-o",dest="out_dir",required=False, type=str, default="results/")
args = vars(parser.parse_args())

bands_file = args["file"]
bands_index = args["band"]
iso_energy = args["iso_energy"]
out_dir = args["out_dir"]

name_field_band = 'band_' + str(bands_index)


# create a new 'CSV Reader'
bZ_BANDS_SIcsv = CSVReader(registrationName='BZ_BANDS_SI.csv', FileName=[
                           bands_file])

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
bZ_BANDS_SIcsvDisplay = Show(
    bZ_BANDS_SIcsv, spreadSheetView1, 'SpreadSheetRepresentation')

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=None, hint=0)

# Properties modified on bZ_BANDS_SIcsvDisplay
bZ_BANDS_SIcsvDisplay.Assembly = ''

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(
    registrationName='TableToPoints1', Input=bZ_BANDS_SIcsv)
tableToPoints1.XColumn = name_field_band
tableToPoints1.YColumn = name_field_band
tableToPoints1.ZColumn = name_field_band

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'kx'
tableToPoints1.YColumn = 'ky'
tableToPoints1.ZColumn = 'kz'

# show data in view
tableToPoints1Display = Show(
    tableToPoints1, spreadSheetView1, 'SpreadSheetRepresentation')

# hide data in view
Hide(bZ_BANDS_SIcsv, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# Properties modified on tableToPoints1Display
tableToPoints1Display.Assembly = ''

# create new layout object 'Layout #2'
layout2 = CreateLayout(name='Layout #2')

# set active view
SetActiveView(None)

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

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView1_1, layout=layout2, hint=0)

# set active source
SetActiveSource(tableToPoints1)

# show data in view
tableToPoints1Display_1 = Show(
    tableToPoints1, renderView1_1, 'GeometryRepresentation')

# trace defaults for the display properties.
tableToPoints1Display_1.Representation = 'Surface'
tableToPoints1Display_1.ColorArrayName = [None, '']
tableToPoints1Display_1.SelectTCoordArray = 'None'
tableToPoints1Display_1.SelectNormalArray = 'None'
tableToPoints1Display_1.SelectTangentArray = 'None'
tableToPoints1Display_1.OSPRayScaleArray = name_field_band
tableToPoints1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display_1.SelectOrientationVectors = 'None'
tableToPoints1Display_1.ScaleFactor = 0.2
tableToPoints1Display_1.SelectScaleArray = 'None'
tableToPoints1Display_1.GlyphType = 'Arrow'
tableToPoints1Display_1.GlyphTableIndexArray = 'None'
tableToPoints1Display_1.GaussianRadius = 0.01
tableToPoints1Display_1.SetScaleArray = ['POINTS', name_field_band]
tableToPoints1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.OpacityArray = ['POINTS', name_field_band]
tableToPoints1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToPoints1Display_1.ScaleTransferFunction.Points = [
    -12.5891, 0.0, 0.5, 0.0, -8.1794, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToPoints1Display_1.OpacityTransferFunction.Points = [
    -12.5891, 0.0, 0.5, 0.0, -8.1794, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1_1.ResetCamera(False)

# reset view to fit data
renderView1_1.ResetCamera(False)


# set scalar coloring
ColorBy(tableToPoints1Display_1, ('POINTS', name_field_band))

# rescale color and/or opacity maps used to include current data range
tableToPoints1Display_1.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
tableToPoints1Display_1.SetScalarBarVisibility(renderView1_1, True)

# get color transfer function/color map for name_field_band
band_0LUT = GetColorTransferFunction(name_field_band)

# get opacity transfer function/opacity map for name_field_band
band_0PWF = GetOpacityTransferFunction(name_field_band)

# reset view to fit data bounds
renderView1_1.ResetCamera(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, False)

SaveScreenshot('k_points_bz_mesh.png', renderView1_1, ImageResolution=[1900, 1500],
               FontScaling="Do not scale fonts",
               OverrideColorPalette='PrintBackground',
               TransparentBackground=0,
               # PNG options
               CompressionLevel='2')

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(registrationName='Delaunay3D1', Input=tableToPoints1)

# show data in view
delaunay3D1Display = Show(delaunay3D1, renderView1_1,
                          'UnstructuredGridRepresentation')

# trace defaults for the display properties.
delaunay3D1Display.Representation = 'Surface'
delaunay3D1Display.ColorArrayName = ['POINTS', name_field_band]
delaunay3D1Display.LookupTable = band_0LUT
delaunay3D1Display.SelectTCoordArray = 'None'
delaunay3D1Display.SelectNormalArray = 'None'
delaunay3D1Display.SelectTangentArray = 'None'
delaunay3D1Display.OSPRayScaleArray = name_field_band
delaunay3D1Display.OSPRayScaleFunction = 'PiecewiseFunction'
delaunay3D1Display.SelectOrientationVectors = 'None'
delaunay3D1Display.ScaleFactor = 0.2
delaunay3D1Display.SelectScaleArray = 'None'
delaunay3D1Display.GlyphType = 'Arrow'
delaunay3D1Display.GlyphTableIndexArray = 'None'
delaunay3D1Display.GaussianRadius = 0.01
delaunay3D1Display.SetScaleArray = ['POINTS', name_field_band]
delaunay3D1Display.ScaleTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.OpacityArray = ['POINTS', name_field_band]
delaunay3D1Display.OpacityTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.DataAxesGrid = 'GridAxesRepresentation'
delaunay3D1Display.PolarAxes = 'PolarAxesRepresentation'
delaunay3D1Display.ScalarOpacityFunction = band_0PWF
delaunay3D1Display.ScalarOpacityUnitDistance = 0.03884480555551276
delaunay3D1Display.OpacityArrayName = ['POINTS', name_field_band]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
delaunay3D1Display.ScaleTransferFunction.Points = [
    -12.5891, 0.0, 0.5, 0.0, -8.1794, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
delaunay3D1Display.OpacityTransferFunction.Points = [
    -12.5891, 0.0, 0.5, 0.0, -8.1794, 1.0, 0.5, 0.0]

# hide data in view
Hide(tableToPoints1, renderView1_1)

# show color bar/color legend
delaunay3D1Display.SetScalarBarVisibility(renderView1_1, True)

# update the view to ensure updated data information
renderView1_1.Update()

SaveScreenshot('k_delaunay_bz_mesh.png', renderView1_1, ImageResolution=[1900, 1500],
               FontScaling="Do not scale fonts",
               OverrideColorPalette='PrintBackground',
               TransparentBackground=0,
               # PNG options
               CompressionLevel='2')

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=delaunay3D1)
contour1.ContourBy = ['POINTS', name_field_band]
contour1.Isosurfaces = [iso_energy]
contour1.PointMergeMethod = 'Uniform Binning'

# find source
tableToPoints1 = FindSource('TableToPoints1')

# find source
bZ_BANDS_SIcsv = FindSource('BZ_BANDS_SI.csv')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for name_field_band
band_0LUT = GetColorTransferFunction(name_field_band)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', name_field_band]
contour1Display.LookupTable = band_0LUT
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = name_field_band
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.1421229572214447
contour1Display.SelectScaleArray = name_field_band
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = name_field_band
contour1Display.GaussianRadius = 0.0071061478610722345
contour1Display.SetScaleArray = ['POINTS', name_field_band]
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', name_field_band]
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [iso_energy]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [iso_energy]

# rescale color and/or opacity maps used to exactly fit the current data range
contour1Display.RescaleTransferFunctionToDataRange(False, True)

# get color transfer function/color map for 'band_0'
band_0LUT = GetColorTransferFunction(name_field_band)

# get opacity transfer function/opacity map for 'band_0'
band_0PWF = GetOpacityTransferFunction(name_field_band)

# hide data in view
# Hide(delaunay3D1, renderView1)
# change representation type
delaunay3D1Display.SetRepresentationType('Feature Edges')

renderView1_1.ResetCamera(True)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
contour1Display.SetScalarBarVisibility(renderView1_1, False)
SaveScreenshot(out_dir + '/' + name_field_band + '_iso_surface' + str(iso_energy) + 'eV_band.png', renderView1_1, ImageResolution=[1900, 1500],
               FontScaling="Do not scale fonts",
               OverrideColorPalette='PrintBackground',
               TransparentBackground=0,
               # PNG options
               CompressionLevel='2')

# get opacity transfer function/opacity map for name_field_band
band_0PWF = GetOpacityTransferFunction(name_field_band)

# # create a new 'Contour'
# contour1 = Contour(registrationName='Contour1', Input=delaunay3D1)
# name_field_band = 'band_' + str(bands_index)
# contour1.ContourBy = ['POINTS', name_field_band]
# contour1.Isosurfaces = [-10.38425]
# contour1.PointMergeMethod = 'Uniform Binning'

# # Properties modified on contour1
# contour1.ContourBy = ['POINTS', name_field_band]
# contour1.Isosurfaces = [-10.38425, 1.03604, 1.465957777777778, 1.8958755555555558, 2.3257933333333334,
#                         2.7557111111111112, 3.185628888888889, 3.615546666666667, 4.045464444444445, 4.475382222222223, 4.9053]

# # show data in view
# contour1Display = Show(contour1, renderView1_1, 'GeometryRepresentation')

# # get color transfer function/color map for name_field_band
# band_4LUT = GetColorTransferFunction(name_field_band)

# # trace defaults for the display properties.
# contour1Display.Representation = 'Surface'
# contour1Display.ColorArrayName = ['POINTS', name_field_band]
# contour1Display.LookupTable = band_4LUT
# contour1Display.SelectTCoordArray = 'None'
# contour1Display.SelectNormalArray = 'Normals'
# contour1Display.SelectTangentArray = 'None'
# contour1Display.OSPRayScaleArray = name_field_band
# contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# contour1Display.SelectOrientationVectors = 'None'
# contour1Display.ScaleFactor = 0.2
# contour1Display.SelectScaleArray = name_field_band
# contour1Display.GlyphType = 'Arrow'
# contour1Display.GlyphTableIndexArray = name_field_band
# contour1Display.GaussianRadius = 0.01
# contour1Display.SetScaleArray = ['POINTS', name_field_band]
# contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
# contour1Display.OpacityArray = ['POINTS', name_field_band]
# contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
# contour1Display.DataAxesGrid = 'GridAxesRepresentation'
# contour1Display.PolarAxes = 'PolarAxesRepresentation'

# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# contour1Display.ScaleTransferFunction.Points = [
#     1.465957777777778, 0.0, 0.5, 0.0, 4.9053, 1.0, 0.5, 0.0]

# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# contour1Display.OpacityTransferFunction.Points = [
#     1.465957777777778, 0.0, 0.5, 0.0, 4.9053, 1.0, 0.5, 0.0]

# # hide data in view
# Hide(delaunay3D1, renderView1_1)

# # show color bar/color legend
# contour1Display.SetScalarBarVisibility(renderView1_1, True)

# # update the view to ensure updated data information
# renderView1_1.Update()

# # get opacity transfer function/opacity map for name_field_band
# band_4PWF = GetOpacityTransferFunction(name_field_band)

# ================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
# ================================================================

# get layout
layout1 = GetLayoutByName("Layout #1")

# --------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(400, 400)

# layout/tab size in pixels
layout2.SetSize(1419, 793)

# -----------------------------------
# saving camera placements for views

# current camera placement for renderView1_1
renderView1_1.CameraPosition = [-1.219498509516459,
                                0.10737148019145519, 4.403821714274042]
renderView1_1.CameraViewUp = [0.003207918728300859,
                              0.9997190131803833, -0.023486250085388804]
renderView1_1.CameraParallelScale = 1.7320508075688772


# save screenshot



# --------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
