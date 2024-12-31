# Overview
A collection of tools in KLayout for navigation, polygon positioning, and photonic waveguide routing

# Navigation/Selection

Navigation shortcuts not available in base KLayout including

- Change drawing/snapping grid, equivalent to adjusting Editor Options>Snapping>Grid (Snapper>Editing>Grid)
- Toggle Top Level Selection (Snapper>Selection>Edit Top Level Selection On/Off)
- Select instances in current layout from cell hierarchy (Snapper>Selection>Select instances)
- Ascent/descent cell instance as new top cell (Snapper>Navigation>Enter Instance, Exit Instance, To Top Cell)
- ...

# Snap

Snapper>Snap>Snap object aligns the first selected polygons/pcells/boxes to the last selection. If there is only 2 selections, the first selected object is aligned the the second on a North-East-South-West basis, between the closest vertical/horizontal parallel faces of the objects. If there are 3 selections, the first selected object will be aligned to the second along one axis and the aligned to the other object along the other axis. Three versions of the function are available for center, left, and right justification of the alignment.

Snapper>Snap>Snap path aligns a selected path to the second selected object, similar to Snap object.

# Waveguide to Path
Snapper>Editing>Routing>Path to waveguide requires a path to be drawn and selected, and choice of a bend cell to replace corners selected through 
Snapper Options (Snapper>Snapper Options menu) in the Waveguide Defaults section.
![image](https://github.com/user-attachments/assets/cfa41a26-5742-4c07-acbb-65530878824a)


Layer, Bend Cell, Width and Radius needs to be configured according to the target replacement bend cell. The origin of the bend cell should be on the bottom left, with respect to its centerline. This function only supports a single layer, single width waveguide at present.

The Path to waveguide function will convert a selected path to a collection of straight segments and bend cells. Bends that do not have sufficient clearance will remain as Path. The original path will remain in place, after the application of the function, to aid reversal of the process if edits are required. A waveguide constructed from Path to waveguide can be cleared via Snapper>Editing>Routing>Revert PCell to path.
