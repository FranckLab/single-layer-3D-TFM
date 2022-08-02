# single-layer-3D-TFM
Repository for single-layer 3D TFM scripts used for measuring cellular displacements and tractions from a gel with a single bead layer imaged with epifluorescence microscopy. This contains the Matlab, Python, and PowerShell scripts needed to run our single-layer 3D TFM algorithm. This algorithm uses epifluorescence image stacks with a single layer of fluorescent micro-beads embedded at a gel surface to reconstuct 3D surface layer displacements and tractions, given approperiate material properties for the gel. See our paper describing this process for more details: Hazlett, L., Landauer, A.K., Patel, M. et al. Epifluorescence-based three-dimensional traction force microscopy. Sci Rep 10, 16599 (2020). https://doi.org/10.1038/s41598-020-72931-6



### Important pages
* [Download latest version!](https://github.com/FranckLab/single-layer-3D-TFM/releases)
* [Example data](https://minds.wisconsin.edu/handle/1793/79908)
* [FAQ](https://github.com/FranckLab/single-layer-3D-TFM#faq)
* [Questions/Issues](https://github.com/FranckLab/single-layer-3D-TFM/issues)
* [Cite](https://github.com/FranckLab/single-layer-3D-TFM#cite)
* [Franck Lab](https://www.franck.engr.wisc.edu)
 
## Running Single Layer 3D TFM

Notes: 
  - For example datasets (both real experiments and sythetic data) see our data storage repository for this project hosted by UW MINDS here: https://minds.wisconsin.edu/handle/1793/79908
  - For more detailed guidelines see the user manual (SI) in our associated paper.

### Input 3D Image Stack Requirements
* To check if the 3D image stack is approperiate, consult the orginal TPT guidelines and our associated [TPT paper](https://github.com/FranckLab/T-PT).
* The 3D image stack should be saved and organized as outlined in the Supplimental Information in our associated paper. 

### Running included example case
0. Install and configure Matlab (with Image Processing and Statistics toolboxes) and Docker Desktop 
1. Make sure that the main files are added to the path in Matlab and the data directory is configured correctly.
2. Download and save the [example data](https://minds.wisconsin.edu/handle/1793/79908) in the example folder. 
3. Run the `SL3DTFM_runfile.m` script in Matlab, then start Docker and run the `run_sl_tfm_ps.ps1` in PowerShell, and finally run the `SL3DTFM_postFEniCS.m` script in Matlab.

### Health warning!
Deconvolution and Docker both may require a **large amount of RAM**. We recommend a minimum of 32GB system RAM. To avoid issues, run the Matlab displacement computations, then start Docker and run the FEA. This works best since at least 16GB of RAM should be allocated to Docker and deconvolution may also require more than 16GB RAM depending on settings and image sizes.

## 3rd party files
We have used several 3rd party Matlab scripts, cheifly from the Matlab File Exchange, and have included relevent licences in a subfolder.
* Supplemental `.m`-files from the MATLAB file exchange include:
 - `gridfit.m`
 - `imagesc3D.m`
 - `inpaint_nans3.m`
 - `regionprops3d.m`
 - `regularizeNd.m`
 - `turbo.m`
 
## Tips and FAQ

* How do I? ...
 - Note that after downloading data from MINDS as a .zip, directly extracting the folder in place will result in an extra folder "layer". Make sure the folder structure is as described in the documentation before running the code. 
 -
 -


## Cite
If used please cite:

```bibtex

@article{hazlett_epifluorescence-based_2020,
	title = {Epifluorescence-based three-dimensional traction force microscopy},
	volume = {10},
	copyright = {2021 The Author(s)},
	issn = {2045-2322},
	url = {https://www.nature.com/articles/s41598-020-72931-6},
	doi = {10.1038/s41598-020-72931-6},
	number = {1},
	journal = {Scientific Reports},
	author = {Hazlett, Lauren and Landauer, Alexander K. and Patel, Mohak and Witt, Hadley A. and Yang, Jin and Reichner, Jonathan S. and Franck, Christian},
	month = oct,
	year = {2020},
	pages = {16599}
}

```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/single-layer-3D-TFM#faq) and [Questions/Issues](https://github.com/FranckLab/single-layer-3D-TFM/issues). Add a new question if similar issue has not yet been reported. The authors' contact information can be found at [Franck Lab](https://www.franck.engr.wisc.edu) or in the linked open-access paper.

