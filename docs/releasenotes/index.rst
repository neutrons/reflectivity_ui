.. release_notes


Release Notes
=============
Notes for major and minor releases. Notes for Patch releases are deferred.


<Next Major or Minor Release>
--------------
(date of release, format YYYY-MM-DD)

**Of interest to the User**:

- MR #XYZ: one-liner description

**Of interest to the Developer:**

- MR #XYZ: one-liner description


v3.2.0
------
2023-11-04

**Of interest to the User**:

- MR #96  drop functionality to save files in genx format
- MR #94  fallback option for loading old data runs
- MR #93  change the default size of the smoothing area, and to change the color map on the preview plot
- MR #77  Clearly distinguish global and per-run options
- MR #75  Export direct beam data
- MR #74  Add plot toolbar buttons for log x-axis, log y-axis and R*Q^4
- MR #73  Make messages in status bar disappear after 10 seconds


**Of interest to the Developer:**

- MR #88  Recalculate and replot reflectivity when global config is changed
- MR #79  Update codecov version and add token in actions.yml
- MR #78  Add ConfigurationHandler which connects signals from global configuration parameter UI widgets to callbacks
- MR #71  disable automatic Conda deployments


v3.1.0
------
2023-10-13

**Of interest to the User**:

- MR #67  Adjust the reduction table to make all data visible
- MR #66  Display error message when "Normalize to unity when stitching" fails
- MR #64  Stitching by fitting a polynomial
- MR #63  Bugfix: Losing direct beam run
- MR #61  Correct error bars of stitched reflectivity curves
- MR #60 Publish documentation
- MR #58 Global stitching
- MR #52 Bugfix error in the x vs TOF plot for empty reflectivity curves
- MR #51 Fix plot of smooth off-specular area
- MR #48 Fix broken IO when saving off-specular data


**Of interest to the Developer:**

- MR #57 Remove class ApplicationConfiguration
- MR #54 Unit testing for stitching
- MR #46 Update enviroment.yml with mantidworkbench


v3.0.0
------
2023-04-25
