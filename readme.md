# Dynamic Kicks (`DYNK` input block)

>*NOTE THAT THE DYNK BLOCK AS DESCRIBED HERE IS NOT YET AVAILABLE IN THE STANDARD SIXTRACK DISTRIBUTION* 

This repository contains the code adding the functionality of dynamic kicks to [SixTrack](http://sixtrack.web.cern.ch/SixTrack/).

The `DYNK` block allows the user to modify the original attributes of the elements of the lattice of the accelerator, per turn, with a user defined function.

* The *function* is defined by the user in the _fort.3_ input file, inside the `DYNK` block, via the `FUN` flag.
* The *element's attribute* to be changed (and the turns when this is applicable) is defined by the user in the _fort.3_ input file, inside the `DYNK` block, via the `SET` flag.

Multiple `FUN` functions can be defined. They may also depend on the functions defined above them in the `DYNK` block. Similarly, multiple `SET` may be defined for multiple element/attribute combinations, and for the same element/attribute given that they are not active on the same turns.

For more information about the syntax and how to use it visit the [SixTrack Twiki](https://twiki.cern.ch/twiki/bin/view/LHCAtHome/SixTrackDoc#Dynamic_Kicks_DYNK_input_block).


### Dependance of subroutines and some functions
* DATEN
	* initialize_element(i,.true.)
	
* TRAUTHIN
	* dynk_pretrack
		* dynk_get_value_single
			* dynk_getvalue
* TRAUTHCK
	* dynk_pretrack
		* dynk_get_value_single
			* dynk_getvalue
* THIN4D
	* dynk_apply
		* dynk_setvalue
			* initialize_element(ii,.false.)
			* dynk_computeFUN
				* dynk_lininterp
		* dynk_getvalue

* THIN6D
	* dynk_apply
		* dynk_setvalue
			* initialize_element(ii,.false.)
			* dynk_computeFUN
				* dynk_lininterp
		* dynk_getvalue

* THIN6DUA
	* dynk_apply
			* dynk_setvalue
			* initialize_element(ii,.false.)
			* dynk_computeFUN
				* dynk_lininterp
		* dynk_getvalue

* THCK4D
	* dynk_apply
			* dynk_setvalue
			* initialize_element(ii,.false.)
			* dynk_computeFUN
				* dynk_lininterp
		* dynk_getvalue

* THCK6D
	* dynk_apply
	 		* dynk_setvalue
			* initialize_element(ii,.false.)
			* dynk_computeFUN
				* dynk_lininterp
		* dynk_getvalue

* THCK6DUA
	* dynk_apply
			* dynk_setvalue
			* initialize_element(ii,.false.)
			* dynk_computeFUN
				* dynk_lininterp
		* dynk_getvalue
