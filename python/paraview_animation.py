# import the simple module from the paraview
from paraview.simple import *
import numpy as np
import argparse
# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


def plot_iso_surface(filename_vtu, band_index, min_energy, max_energy, number_iso_values, out_file, nb_frames):

    list_energies = np.linspace(min_energy, max_energy, number_iso_values)
    band_str = f"band_{band_index}"
    band_color = f"band_{(band_index+1)%12}"

    # create a new 'XML Unstructured Grid Reader'
    medium_bandsvtu = XMLUnstructuredGridReader(registrationName='medium_bands.vtu', FileName=['/home/remi/EmpiricalPseudopotential/build/medium_bands.vtu'])
    medium_bandsvtu.CellArrayStatus = []
    medium_bandsvtu.PointArrayStatus = ['band_0', 'band_1', 'band_2', 'band_3', 'band_4', 'band_5', 'band_6', 'band_7', 'band_8', 'band_9', 'band_10', 'band_11']
    medium_bandsvtu.TimeArray = 'TimeValue'

    # Properties modified on medium_bandsvtu
    medium_bandsvtu.TimeArray = 'None'

    UpdatePipeline(time=0.0, proxy=medium_bandsvtu)

    # set active source
    SetActiveSource(medium_bandsvtu)

    UpdatePipeline(time=0.0, proxy=medium_bandsvtu)

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=medium_bandsvtu)
    contour1.ContourBy = ['POINTS', 'band_0']
    contour1.ComputeNormals = 1
    contour1.ComputeGradients = 0
    contour1.ComputeScalars = 1
    contour1.OutputPointsPrecision = 'Same as input'
    contour1.GenerateTriangles = 1
    contour1.Isosurfaces = [-10.38426347869867]
    contour1.PointMergeMethod = 'Uniform Binning'

    # init the 'Uniform Binning' selected for 'PointMergeMethod'
    contour1.PointMergeMethod.Divisions = [50, 50, 50]
    contour1.PointMergeMethod.Numberofpointsperbucket = 8

    renderView1 = GetActiveViewOrCreate('RenderView')
    UpdatePipeline(time=0.0, proxy=contour1)

    # get animation track
    contour1ContourValuesTrack = GetAnimationTrack('ContourValues', index=-1, proxy=contour1)

    # create keyframes for this animation track

    # create a key frame
    keyFrame9974 = CompositeKeyFrame()
    keyFrame9974.KeyTime = 0.0
    keyFrame9974.KeyValues = [-12.58912229538956]
    keyFrame9974.Interpolation = 'Ramp'
    keyFrame9974.Base = 2.0
    keyFrame9974.StartPower = 0.0
    keyFrame9974.EndPower = 1.0
    keyFrame9974.Phase = 0.0
    keyFrame9974.Frequency = 1.0
    keyFrame9974.Offset = 0.0

    # create a key frame
    keyFrame9975 = CompositeKeyFrame()
    keyFrame9975.KeyTime = 1.0
    keyFrame9975.KeyValues = [-8.17940466200778]
    keyFrame9975.Interpolation = 'Ramp'
    keyFrame9975.Base = 2.0
    keyFrame9975.StartPower = 0.0
    keyFrame9975.EndPower = 1.0
    keyFrame9975.Phase = 0.0
    keyFrame9975.Frequency = 1.0
    keyFrame9975.Offset = 0.0

    # initialize the animation track
    contour1ContourValuesTrack.TimeMode = 'Normalized'
    contour1ContourValuesTrack.StartTime = 0.0
    contour1ContourValuesTrack.EndTime = 1.0
    contour1ContourValuesTrack.Enabled = 1
    contour1ContourValuesTrack.KeyFrames = [keyFrame9974, keyFrame9975]
    
    cameraAnimationCue1 = GetCameraTrack(view=renderView1)
    # initialize the animation track
    cameraAnimationCue1.Mode = 'Path-based'
    cameraAnimationCue1.KeyFrames = [keyFrame9974, keyFrame9975]
    
    
    # save animation
    SaveAnimation("TEST.avi", renderView1, ImageResolution=[1200, 1200],
        FontScaling='Do not scale fonts',
        FrameRate=12,
        Compression=True,
        Quality=2, 
        FrameWindow=[0, 100])




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
