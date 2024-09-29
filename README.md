# HeLab_CalciumWave_Simulation_2024
This is a reposity of matlab codes for the work: Sensing of dietary amino acids and regulation of calcium dynamics in adipose tissues through Adipokinetic hormone in Drosophila. The code was built by Shengyao Luo. 
## Instructions to run the code on a personal machine
1. We strongly recommend you to download the codefile on your personal machine to run. You can download all of the associated code [here](https://github.com/Sy-Luo/HeLab_CalciumWave_Simulation_2024/archive/refs/heads/main.zip).
- Each notebook is labeled corresponding to its associated figure in the main text of the manuscript.
- NOTE: The full folder is approximately 10 MB of data as it contains the pre-run simulation output pictures and videos.
2. Extract the downloaded .zip folder where you would like to run the code.
- You should locate all the file in the same document.
3. The version of MATLAB we used is [MATLAB R2023a-academic use](https://www.mathworks.cn/products/matlab.html). You can use the more updated version to run the codes.
4. Click the App-Ons and make sure that you have downloaded the following toolboxs:
  ![](/readme_fig/rm_fig1.png)
  - Curve Fitting Toolbox (ver. 4.17)
  - Image Aquisition Toolbox (ver. 6.7.1)
  - Image Processing Toolbox (ver. 11.7)
  - MATLAB Coder (ver. 5.6)
  - Optimization Toolbox (ver. 9.5)
  - Parallel Computing Toolbox (ver. 7.8)
  - Statistics and Machine Learning Toolbox (ver. 12.5)
  - Symbolic Math Toolbox (ver. 9.3)
5. Use MATLAB to open each individual .m file to run the simulations as desired.

## Code for repeated simulations
You can find the code and corresponding outputs for repeated simulations (with different random number generator seeds) for the manuscript's main figures in the following locations:
- [random_wave_property.m](/random_wave_property.m)

## Result files
You can find the .avi files as the representative results of the codes of the same name:
- [adult_ICW_wildtype.avi](/adult_ICW_wildtype.avi)
- [adult_ICW_inx2i.avi](/adult_ICW_inx2i.avi)
- [larva_ICW_wildtype.avi](/larva_ICW_wildtype.avi)
- [larva_ICW_inx2i.avi](/larva_ICW_inx2i.avi)

## Repository last updated: September 30, 2024 0:00AM GMT+8
