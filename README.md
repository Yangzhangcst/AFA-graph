# AFA-Graph
 
The official codes for "Adaptive Fusion Affinity Graph with Noise-free Online Low-rank Representation for Natural Image Segmentation"
The AFA-Graph paper can be found [here](https://arxiv.org/abs/2110.11685). A list of papers and datasets about natural/color image segmentation can be found in our [repository](https://github.com/Yangzhangcst/Natural-color-image-segmentation).

AFA-Graph is modified from [the offcial SAS implementation](http://www.ee.columbia.edu/ln/dvmm/SuperPixelSeg/dlform.htm).


### Requirements
The code requires the version of Matlab2018a, Ubuntu 16.04.

### Data
The original MSRC dataset can be downloaded from [here](https://www.microsoft.com/en-us/research/project/image-understanding/?from=http%3A%2F%2Fresearch.microsoft.com%2Fvision%2Fcambridge%2Frecognition%2F). We place the processed dataset in `database/MSRC/` folder.

### Demo
Run the demo `IKDE_BSDS300_searchN.m` for 'BSD300'.  
Run the demo `IKDE_BSDS500_searchN.m` for 'BSD500'.  
Run the demo `IKDE_MSRC_searchN.m` for 'MSRC'.  
Run the demo `IKDE_SBD_searchN.m` for 'SBD'.  
Run the demo `IKDE_VOC_searchN.m` for 'PASCAL VOC'.

### Results of FIVE datasets in the paper
The detailed results of BSD300 can be found in `results_searchN_denoise_n/ikde_1_1/evaluation_BSD300.txt`.  
The detailed results of BSD500 can be found in `results_searchN_denoise_n/ikde_1_1/evaluation_BSD500.txt`.  
The detailed results of MSRC can be found in `results_searchN_denoise_n/ikde_1_1/evaluation_MSRC.txt`.  
The detailed results of SBD can be found in `results_searchN_denoise_n/ikde_1_1/evaluation_SBD.txt`.  
The detailed results of VOC can be found in `results_searchN_denoise_n/ikde_1_1/evaluation_VOC.txt`.

### Citing
If you find this repository useful in your research, please consider citing:
```
@INPROCEEDINGS{AFA-Graph,  
  author={Y. {Zhang} and M. {Liu} and H. {Zhang} and G. {Sun} and J. {He}},  
  booktitle={arXiv:2110.11685},   
  title={Adaptive Fusion Affinity Graph with Noise-free Online Low-rank Representation for Natural Image Segmentation},   
  year={2021}}
```
