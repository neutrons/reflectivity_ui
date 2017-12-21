# Magnetic Reflectivity Reduction

# Interface
- Make sure matplotlib widgets work.
- Add progress bar
- Add log window
- Add pop-up message boxes for errors?
- Create a widget to change the mantid path
- Add Huber X cut
- Update affected reflectivity when a direct beam run is deleted or changed.

## Data Manager Design
- Improve data manager design to handle direct beams more cleanly
- Keep the Mantid workspaces instead of deleting them so we can recalculate faster.
- Figure out which cross-section to use for the direct beam if there are more than a single one.
- Keep the Mantid workspaces instead of deleting them so we can recalculate faster.
- Deal with peak position when it's not in the middle of the range
- Add filter to direct beam data too.
- Allow for direct beam data to be a workspace instead of just a run number.

## QUESTION: should we use the same ranges for all cross-sections? Probably.

## Reduction
- Stitching
- Normalize output to 1.

## DONE:
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


