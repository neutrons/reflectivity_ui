# Magnetic Reflectivity Reduction
To-do list:

## Interface
- Make sure matplotlib widgets work.
- Add progress bar
- Add log window
- Add pop-up message boxes for errors?
- Create a widget to change the mantid path
- Add Huber X cut
- Update affected reflectivity when a direct beam run is deleted or changed.
- Do something with DIRPIX and DANGLE0 overwrite
- Make sure all the config properties make it to QSettings
- Refactor the configuration object of the NexusData class
- Add mantid script as a possible output file
- Jupyter notebook as output
- Proper mantid script as output
- Use PolarizerLabel and AnalyzerLabel to better decide on cross-section names
- Use Analyzer log to decide if the analyzer was used and name cross-sections accordingly.
- Plotting: allow for log(x)
- Plot slices in 2D off-specular data

## Data Manager Design
- Improve data manager design to handle direct beams more cleanly
- Figure out which cross-section to use for the direct beam if there are more than a single one.
- Deal with peak position when it's not in the middle of the range

## Reduction
- Stitching  -> make Stitch1DMany take scaling factor for the first workspace
- Trim workspaces when stitching
- Add option for constant Q
- Add an option to choose which cross-sections to use for asymmetry calculation
- Add Q binning option for output
- Background subtraction for offspec and GISANS
- Mantid: add option not to rebin in Q at the end.

## Testing
- Test numpy output
- Test Matlab output
- Test Genx output

## DONE
- Normalize output to 1 (broken when turned off).
- Norm TOF-lambda doesn't work.
- Do something with P0 / PN to clean up the reflectivity
- Add indicator for whether the application things the data set we are looking at is a direct beam or not.
- Decide whether to keep peak ranges common to all cross-sections.
- Switch active data when clicking on the reduction table.
- Tie main UI to normalization list as well as reduction list.
- When clicking on an entry in the reduction table:
 - Flag as active with a color indicator.
 - Get it's config and push it to the main UI.

- When loading a data set from cache:
 - Check whether it's also in the reduction table.
 - If so, get it's config and push it to the main UI.

- When updating the main UI, check whether we need to update the reduction table.
- Match direct beam should update the reduction list.
- Implement item handling in the direct beam list to be identical to the scattering run list.
- Calculated sangle needs to update if we change the position ranges
- Keep the Mantid workspaces instead of deleting them so we can recalculate faster.
- Add filter to direct beam data too.
- Allow for direct beam data to be a workspace instead of just a run number.
- Add InputNormalizationWorkspace as input to MagnetismReflectometerReduction
- Fixed saving/loading of peak ranges
- Add option to cut the first N and last M points -> Needs to be saved properly to file


