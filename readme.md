# Dynamic Kicks (`DYNK` input block)

>*NOTE THAT THE DYNK BLOCK AS DESCRIBED HERE IS NOT YET AVAILABLE IN THE STANDARD SIXTRACK DISTRIBUTION* 

This repository contains the code adding the functionality of dynamic kicks to [SixTrack](http://sixtrack.web.cern.ch/SixTrack/).

The `DYNK` block allows the user to modify the original attributes of the elements of the lattice of the accelerator, per turn. 

This is done with the `SET` flag. The attributes will change following the functions defined by the user, using the `FUN` flag. Multiple `FUN` functions can be defined. They may also depend on the functions defined above them in the `DYNK` block. Similarly, multiple `SET` may be defined for multiple element/attribute combinations, and for the same element/attribute given that they are not active on the same turns.

Additionally, giving the `DEBU` flag turns on extra output to stdout, and lines with "/" as the first character are treated as comments and ignored.

For more information visit the [SixTrack Twiki](https://twiki.cern.ch/twiki/bin/view/LHCAtHome/SixTrackDoc#Dynamic_Kicks_DYNK_input_block).
