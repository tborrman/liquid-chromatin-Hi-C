# Liquid-Chromatin-Hi-C
![LC-HiC](https://github.com/tborrman/liquid-chromatin-Hi-C/blob/master/figures/LChic_progression.PNG)

All code developed for the analysis of Liquid-Chromatin Hi-C assays. 

A majority of code performs downstream analysis of Hi-C interaction matrices in cworld compatible HDF5 format:
- https://github.com/dekkerlab/cMapping
- https://github.com/dekkerlab/cworld-dekker

with ongoing devlopment for compatibility with cooler format:
- https://github.com/mirnylab/cooler
## Usage Examples
### Calculating Loss of Structure (LOS)
```bash
# Scale control Hi-C matrix and predigested Hi-C matrix
src/liquid_chromatin_HiC/scale.py -i control.hdf5
src/liquid_chromatin_HiC/scale.py -i predigest.hdf5

# Calculate 6 Mb cis percentage for scaled matrices
src/liquid_chromatin_HiC/hdf2RangeCisPercent.py -i control_scaleBy_x.hdf5
src/liquid_chromatin_HiC/hdf2RangeCisPercent.py -i predigest_scaleBy_x.hdf5

# Calculate LOS from cis percent bedGraphs
src/scripts/LOS.R -i predigest_cispercent.bedGraph -m control_cispercent.bedGraph
```
### Calculating t<sub>1/2</sub>
```bash
# With all *LOS.bedGraphs from multiple timepoints in directory
src/scripts/make_LOS_table.R
# With output LOS table of all timepoints
src/scripts/get_half_life_exponential.R -i LOS_40kb.txt -f True
```
## Contributing to liquid-chromatin-Hi-C
All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome!

## References
Houda Belaghzal*, Tyler Borrman*, Andrew D. Stephens, Denis L. Lafontaine, Sergey V. Venev, Zhiping Weng, John F. Marko, Job Dekker. (2019). [Compartment-dependent chromatin interaction dynamics revealed by liquid chromatin Hi-C](https://www.biorxiv.org/content/10.1101/704957v1). *bioRxiv*

## License 
[GNU General Public License v3.0](LICENSE)


