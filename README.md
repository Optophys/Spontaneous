# Spontaneous
Code for "Spontaneous activity competes with externally evoked responses in sensory cortex", accepted for publication in PNAS and archived at https://www.biorxiv.org/content/10.1101/2020.08.18.256206v2.
The matlab code in Karvat_et_al_2021_script_1.m shows the calculations performed to obtain figure 2, with accompanying raw data from one session.
Other m. files are accessory functions which should be added to the path.
For power calculations, we use the fieldtrip toolbox, which has to be added to the path as well. This free tooblox can be downloaded from https://www.fieldtriptoolbox.org/.
.mat files:
The raw LFP from 16 channels over the whole session was split into 4 files (lfp1-4) so they could be uploaded to GitHub. The code in Karvat_et_al_2021_script_1.m concatenates them.
FR contains the z-normalized population firing-rate.
ratP is a structure with the behavioral outcome and timestamps per trial.
