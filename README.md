# Supplementary code repository for GAP XYZ
Here, we provide a basic implementation of the standardization procedure discussed in GAP XYZ.
The method is applicable to all events measured by a detector aranged in an isometric triangular grid, when the azimuth direction is well known.
While doing this it reduces the phase space approximately by a factor of twelve.


## Basic code description
The script XXX.py takes a list of integral position corresponding to the basis defined in GAP XYZ.
It is supposed that the first item in the list is the central point around which the event should be normalized.
The output are the new integral positions corresponding to the same basis.


## Implementation details
The code is written in Python 3.
We solely use numpy as backend.
However, for the plotting routines also matplotlib is used.
