# Magnetic Reflectivity Reduction

#TODO
- Add indicator for whether the application things the data set we are looking at is a direct beam or not.
- Make sure matplotlib widgets work.
- Implement item handling in the direct beam list to be identical to the scattering run list (add a remove button).
- Decide whether to keep peak ranges common to all cross-sections.
- Figure out which cross-section to use for the direct beam if there are more than a single one.
- Add progress bar
- Add log window
- Improve data manager design to handle direct beams more cleanly
- Keep the Mantid workspaces instead of deleting them so we can recalculate faster.
- Switch active data when clicking on the reduction table.
- Calculated sangle needs to update if we change the position ranges