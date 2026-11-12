# Two Photon Data Ana
1. Run plot cnmfe_dff_use.m
   choose neuron.mat(neuron activities analysed by cnmf) & data_64_64_13.mat(use for calculate background)
   ouput&save neuron_dff.mat which means df/f
2. Run export_AccuSleep_State_2.m
   choose neuron_dff.mat
          EEG.mat
          EMG.mat
          labels.mat
   output figure Brain State EEG Spectrum Neuronal Activists
   <img width="2000" height="725" alt="image" src="https://github.com/user-attachments/assets/8a0dff71-503c-418e-82a6-a5f725389fb6" />

## cnmfe_dff_use.m
Note!!!
   Fs?

## export_AccuSleep_State_2.m
Note!!!
   epoch?
