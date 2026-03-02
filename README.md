# cse284-hmm

Implementation of HMM (Local Ancestry Inferene)

Note: if you are using Windows, install WSL so that you can run Linux.
- `wsl --install`
- You may have to restart your computer. If so, after restarting, a terminal window will open automatically to finish the setup.
- You will be asked to create a username and password.
- After setup, you will automatically enter into WSL. To manually enter in the future, use `wsl -d Ubuntu`. To exit WSL, run `exit`.
- Install miniconda in WSL: 
   1. Download the installer script
      `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
   2. Run the installer
      `bash Miniconda3-latest-Linux-x86_64.sh`
   3. Activate bash
      `source ~/.bashrc`

Create an environment:
- In your root directory, run `conda env create -f environment.yml`
= To activate the environment: run `conda activate hmm_env`

## Step 1: Get the Data

1000 Genomes Project Phase 3

1. **Genetic Map**
   - Purpose: to consider the actual biological distances within the chromosome.
2. **Sample Panel**
   - Purpose: contains samples for both populations.
3. **VCF Slice**
   - The part of chromosome 22 that I'm going to test.

To download data:
- `cd data` & `bash ./download_data.sh`
TODO:

- Write script to download the 3 files above (_did not push to git because files are too big_)
- Check sample file

