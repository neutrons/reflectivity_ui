# Magnetic Reflectivity Reduction
To-do list:

## Interface
- Add log window
- Add pop-up message boxes for errors?
- Create a widget to change the mantid path
- Refactor the configuration object of the NexusData class
- Jupyter notebook as output
- Proper mantid script as output
- Plotting: allow for log(x)

## Data Manager Design
- Add option to match direct beam cross-seciont to data cross-section, otherwise sum up all cross-section in the direct beam data file.
- Emptying cache should also delete the mantid workspaces

## Reduction
- Add use-sangle option to output QuickNXS file.
- Stitching  -> make Stitch1DMany take scaling factor for the first workspace
- Add an option to choose which cross-sections to use for asymmetry calculation
- Add Q binning option for output
- Mantid: add option not to rebin in Q at the end.
- Write output algo that takes group workspace.
- Automate selection of const-Q binning

## Off-spec
- Only recalculate off-spec when needed (need to add a need_recalc data member)
- Check that the overall normalization is the same as with old QuickNXS when validating
- Off-spec needs to use option to skip first and last points
  (which only works on specular data, which uses a different binning)
  We may want to harmonize the number of skipped points according to the difference between
  the specular and off-specular binning.
- Try using events directly when computing offspecular.
- Plot slices in 2D off-specular data

## GISANS
- Add background subtraction
- Option to skip first and last points
- Generate binned output

## Testing
- Automated testing
- Test numpy output
- Test Matlab output
- Test Genx output

## DONE
- Add off-spec coordinates to toolbar
- Verify off-spec normalization vs peak ranges and two-theta
- Use PolarizerLabel and AnalyzerLabel to better decide on cross-section names
- Use Analyzer log to decide if the analyzer was used and name cross-sections accordingly.
- Update affected reflectivity when a direct beam run is deleted or changed.
- Make sure all the config properties make it to QSettings
- Add mantid script as a possible output file
- Background subtraction for offspec
- Trim workspaces when stitching
- Add option for constant Q
- Do something with DIRPIX and DANGLE0 overwrite
- Make sure matplotlib widgets work.
- Improve data manager design to handle direct beams more cleanly
- SORT RUNS WHEN LOADING REDUCED FILE
- Add progress bar
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
- Deal with peak position when it's not in the middle of the range

