*  Version 1.5.0:
    * Major release. 
        + Restructure the code: removed Dendropy folder and many redundant files. 
        + Activate installation by Anaconda. 
        + Allow backward time specification via -b flag
        + Perform automatic label assignment and allow suppression by -k flag
        + Allow each calibration point to be specified as LCA of a list of species
*  Version 1.4.0:
    * Major release. Allow calibrations for internal nodes and allow missing sampling times for leaves, hard constraints only. Requires a unique label for each node; sampling times / calibration points are identified by node labels.
*  Version 1.3.1:
    * Added zero-length and verbose options. Remove options for LSD and LF. Remove options for scaling. Remove cvxpy requirement.  
*  Version 1.3.0:
	* Pseudo-count in wLogDate weighting 
* Version 1.2.0:
	* Apply square-root scaling strategy
* Version 1.1.0:
	* Add -e option to remove short terminal branches. Add option to run LF and LSD objective. Use cvxpy with Mosek for LF
* Version 1.0.0:
	* The version submitted to RECOMB 2019.
