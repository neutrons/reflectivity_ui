## v2.0.43 [06/11/2020]
 - Pinned to use mantid42 and python2 fixed
 
## v2.0.42 [05/28/2020]
 - Pinned to use mantid42 and python2

## v2.0.41 [06/29/2019]
 - Minor fixes

## v2.0.40 [02/19/2019]
 - Write calculated angle to data file

## v2.0.39 [02/14/2019]
 - Order cross-sections starting with off-off
 - Fix issue with 3rd cross-section not displaying when only three are measured.

## v2.0.38 [02/12/2019]
 - Fix issue with nan's appearing with fine binning of off-specular

## v2.0.37 [02/11/2019]
 - Improve off-specular slices

## v2.0.36 [02/07/2019]
 - Result preview improvements (sync zoom)
 - Off-spec UI improvements: separate smoothing from the rest
 - Add smoothing dialog for pick options interactively on a plot
 - Email option [data only]
 - When stitching, pick the cross-section with highest counts in overlap
 - Add GISANS result tiles, with sync zoom
 - Add smart stitch that picks the cross-section with highest stats at the overlap
 - Tweak size for smaller screens
 - Display Qz bin width when binning off-spec
 - Improve output python script

## v2.0.35 [01/24/2019]
 - Add result preview

## v2.0.34
 - Make sure spin asymmetry work when not rebinning.
 - Remove the option to use the data's TOF range because it creates problems when combining data (like spin asymmetry).

## v2.0.33
 - Hybrid spin boxes to test out synching to plots and responsiveness.
 - Background selection is now remembered from run to run upon loading.
 - Selecting a very wide wavelength band will pick the actual band in the data.
 - Made sure cutting TOF bins matches between specular (not rebinned) and off-specular
 - Updated wavelength band.
