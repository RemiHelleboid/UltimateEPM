# import the simple module from the paraview
from paraview.simple import *
import numpy as np
import argparse
import inspect
# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


def plot_iso_surface(filename_vtu, band_index, min_energy, max_energy, number_iso_values, out_file, nb_frames):

    list_energies = np.linspace(min_energy, max_energy, number_iso_values)
    band_str = f"band_{band_index}"
    band_color = f"band_{(band_index+1)%12}"

    # create a new 'XML Unstructured Grid Reader'
    medium_1_bz_meshvtu = XMLUnstructuredGridReader(registrationName='fine_1_bz_mesh.vtu', FileName=['/home/remi/EmpiricalPseudopotential/build/fine_1_bz_mesh.vtu'])
    medium_1_bz_meshvtu.PointArrayStatus = ['band_0', 'band_1', 'band_2', 'band_3', 'band_4', 'band_5', 'band_6', 'band_7', 'band_8', 'band_9', 'band_10', 'band_11', 'band_12', 'band_13', 'band_14', 'band_15']


    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # get display properties
    medium_1_bz_meshvtuDisplay = GetDisplayProperties(medium_1_bz_meshvtu, view=renderView1)

    # set scalar coloring
    ColorBy(medium_1_bz_meshvtuDisplay, ('POINTS', 'band_0'))

    # rescale color and/or opacity maps used to include current data range
    medium_1_bz_meshvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    medium_1_bz_meshvtuDisplay.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'band_0'
    band_0LUT = GetColorTransferFunction('band_0')

    # get opacity transfer function/opacity map for 'band_0'
    band_0PWF = GetOpacityTransferFunction('band_0')

    # rescale color and/or opacity maps used to exactly fit the current data range
    medium_1_bz_meshvtuDisplay.RescaleTransferFunctionToDataRange(False, True)

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=medium_1_bz_meshvtu)
    contour1.ContourBy = ['POINTS', 'band_0']
    contour1.Isosurfaces = [-10.379030287270895]
    contour1.PointMergeMethod = 'Uniform Binning'

    # get animation track
    contour1ContourValuesTrack = GetAnimationTrack('ContourValues', index=-1, proxy=contour1)

    # create keyframes for this animation track

    # create a key frame
    keyFrame23966 = CompositeKeyFrame()
    keyFrame23966.KeyValues = [-12.58304167633501]

    # create a key frame
    keyFrame23967 = CompositeKeyFrame()
    keyFrame23967.KeyTime = 1.0
    keyFrame23967.KeyValues = [-8.17501889820678]

    # initialize the animation track
    contour1ContourValuesTrack.KeyFrames = [keyFrame23966, keyFrame23967]

    # get animation scene
    animationScene1 = GetAnimationScene()

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = ['POINTS', 'band_0']
    contour1Display.LookupTable = band_0LUT
    contour1Display.SelectTCoordArray = 'None'
    contour1Display.SelectNormalArray = 'None'
    contour1Display.SelectTangentArray = 'None'
    contour1Display.OSPRayScaleArray = 'band_0'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 1.1564823173178673e-19
    contour1Display.SelectScaleArray = 'band_0'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'band_0'
    contour1Display.GaussianRadius = 5.782411586589336e-21
    contour1Display.SetScaleArray = ['POINTS', 'band_0']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'band_0']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [-8.17501889820678, 0.0, 0.5, 0.0, -8.173066139221191, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [-8.17501889820678, 0.0, 0.5, 0.0, -8.173066139221191, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(medium_1_bz_meshvtu, renderView1)

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    # Properties modified on animationScene1
    animationScene1.NumberOfFrames = 80

    # update the view to ensure updated data information
    renderView1.Update()

    animationScene1.Play()

    # get animation scene
    animationScene1 = GetAnimationScene()

    # save animation
    SaveAnimation('/home/remi/EmpiricalPseudopotential/build/New Folder/bbbe.avi', renderView1, ImageResolution=[1200, 1200],
        FontScaling='Do not scale fonts',
        FrameRate=45,
        FrameWindow=[0, 79])



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', description='Test')
    parser.add_argument("-f", "--filename",dest="file", required=True, type=str)
    parser.add_argument("-b", "--index-band",dest="band", required=False, type=int, default=0)
    parser.add_argument("-m", "--min-energy",dest="min_iso_energy",
                    required=False, type=float, default=-10)
    parser.add_argument("-M", "--max-energy",dest="max_iso_energy",
                    required=False, type=float, default=10)
    parser.add_argument("-n", "--number_iso_values",dest="number_iso_values",
                    required=False, type=int, default=10)
    parser.add_argument("-F", "--nb-frame",dest="nb_frames",
                    required=False, type=int, default=120)
    parser.add_argument("-o", dest="out_dir", required=False,
                    type=str, default="results/")
    args = vars(parser.parse_args())

    bands_file = args["file"]
    bands_index = args["band"]
    min_iso_energy = args["min_iso_energy"]
    max_iso_energy = args["max_iso_energy"]
    number_iso_values = args["number_iso_values"]
    nb_frames = args["nb_frames"]
    out_dir = args["out_dir"]
    
    out_file = f"{out_dir}/rotation_animation_{bands_index}th_band_iso.avi"
    
    plot_iso_surface(bands_file, bands_index, min_iso_energy, max_iso_energy, number_iso_values, out_file, nb_frames)
