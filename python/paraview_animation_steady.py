# import the simple module from the paraview
from paraview.simple import *
import numpy as np
import argparse
import inspect
# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

si_extremas = [[-12.5891, -8.33024],
               [-8.30847, -1.66159e-13],
               [-4.42516, -2.71896e-14],
               [-3.02394, 0.0],
               [1.03669, 3.40913],
               [1.16998, 4.97943],
               [3.37964, 12.2698],
               [3.99022, 12.2698],
               [6.75105, 12.8557],
               [7.31399, 13.7759],
               [8.25655, 16.2184],
               [11.4349, 16.6692],
               [12.569, 17.631],
               [12.569, 18.665],
               [15.5224, 21.5633],
               [17.5607, 27.3715],
               [17.5607, 27.3715],
               [18.4241, 27.421],
               [21.258, 28.5222],
               [21.4927, 29.9422],
               [22.2918, 30.3369],
               [22.2918, 30.9994],
               [26.1788, 33.2535],
               [27.1347, 34.9097],
               [30.3366, 35.8184],
               [31.2036, 36.2169],
               [31.6554, 37.3205],
               [32.5713, 43.4564],
               [34.6604, 43.9431],
               [35.6336, 43.9431]]


def plot_iso_surface(filename_vtu, band_index, min_energy, max_energy, number_iso_values, out_file, nb_frames):

    list_energies = np.linspace(min_energy, max_energy, number_iso_values)
    band_str = f"band_{band_index}"
    band_color = f"band_{(band_index+1)%12}"
    
    min_energy, max_energy = si_extremas[band_index]
    
    print("band_index: ", band_index)
    print("min_energy: ", min_energy)
    print("max_energy: ", max_energy)

    # create a new 'XML Unstructured Grid Reader'
    medium_1_bz_meshvtu = XMLUnstructuredGridReader(registrationName='fine_1_bz_mesh.vtu', FileName=[
                                                    '/home/remi/EmpiricalPseudopotential/build/fine_1_bz_mesh.vtu'])
    medium_1_bz_meshvtu.PointArrayStatus = ['band_0', 'band_1', 'band_2', 'band_3', 'band_4', 'band_5',
                                            'band_6', 'band_7', 'band_8', 'band_9', 'band_10', 'band_11', 'band_12', 'band_13', 'band_14', 'band_15']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # get display properties
    medium_1_bz_meshvtuDisplay = GetDisplayProperties(
        medium_1_bz_meshvtu, view=renderView1)

    # set scalar coloring
    ColorBy(medium_1_bz_meshvtuDisplay, ('POINTS', band_str))

    # rescale color and/or opacity maps used to include current data range
    medium_1_bz_meshvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    medium_1_bz_meshvtuDisplay.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for band_str
    band_0LUT = GetColorTransferFunction(band_str)

    # get opacity transfer function/opacity map for band_str
    band_0PWF = GetOpacityTransferFunction(band_str)

    # rescale color and/or opacity maps used to exactly fit the current data range
    medium_1_bz_meshvtuDisplay.RescaleTransferFunctionToDataRange(False, True)

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=medium_1_bz_meshvtu)
    contour1.ContourBy = ['POINTS', band_str]
    contour1.Isosurfaces = [-10.379030287270895]
    contour1.PointMergeMethod = 'Uniform Binning'

    # Hide orientation axes
    renderView1.OrientationAxesVisibility = 0

    # get animation track
    contour1ContourValuesTrack = GetAnimationTrack(
        'ContourValues', index=-1, proxy=contour1)

    # create keyframes for this animation track

    # create a key frame
    keyFrame23966 = CompositeKeyFrame()
    keyFrame23966.KeyValues = [min_energy]

    # create a key frame
    keyFrame23967 = CompositeKeyFrame()
    keyFrame23967.KeyTime = 1.0
    keyFrame23967.KeyValues = [max_energy]

    # initialize the animation track
    contour1ContourValuesTrack.KeyFrames = [keyFrame23966, keyFrame23967]

    # get animation scene
    animationScene1 = GetAnimationScene()

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = ['POINTS', band_str]
    contour1Display.LookupTable = band_0LUT
    contour1Display.SelectTCoordArray = 'None'
    contour1Display.SelectNormalArray = 'None'
    contour1Display.SelectTangentArray = 'None'
    contour1Display.OSPRayScaleArray = band_str
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 1.1564823173178673e-19
    contour1Display.SelectScaleArray = band_str
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = band_str
    contour1Display.GaussianRadius = 5.782411586589336e-21
    contour1Display.SetScaleArray = ['POINTS', band_str]
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', band_str]
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [
        -8.17501889820678, 0.0, 0.5, 0.0, -8.173066139221191, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [
        -8.17501889820678, 0.0, 0.5, 0.0, -8.173066139221191, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(medium_1_bz_meshvtu, renderView1)

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, False)

    # Properties modified on animationScene1
    animationScene1.NumberOfFrames = nb_frames

    # update the view to ensure updated data information
    renderView1.Update()

    animationScene1.Play()

    # get animation scene
    animationScene1 = GetAnimationScene()

    f
    # save animation
    SaveAnimation(out_file, renderView1, ImageResolution=[1200, 1200],
                  FontScaling='Do not scale fonts',
                  FrameRate=30,
                  FrameWindow=[0, nb_frames-1])


def print_agrs(band_file, min_energy, max_energy, band_index):
    print("band_file: ", band_file)
    print("band_index: ", band_index)
    print("min_energy: ", min_energy)
    print("max_energy: ", max_energy)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', description='Test')
    parser.add_argument("-f", "--filename", dest="file",
                        required=True, type=str)
    parser.add_argument("-b", "--index-band", dest="band",
                        required=False, type=int, default=0)
    parser.add_argument("-m", "--min-energy", dest="min_iso_energy",
                        required=False, type=float, default=-10)
    parser.add_argument("-M", "--max-energy", dest="max_iso_energy",
                        required=False, type=float, default=10)
    parser.add_argument("-n", "--number_iso_values", dest="number_iso_values",
                        required=False, type=int, default=10)
    parser.add_argument("-F", "--nb-frame", dest="nb_frames",
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

    out_file = f"{out_dir}/Si_steady_animation_{bands_index}th_band_iso.avi"

    print_agrs(bands_file, min_iso_energy, max_iso_energy, bands_index)

    plot_iso_surface(bands_file, bands_index, min_iso_energy,
                     max_iso_energy, number_iso_values, out_file, nb_frames)
