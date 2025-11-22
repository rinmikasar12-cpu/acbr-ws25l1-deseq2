# Create Conda environment with Miniconda
I am creating a conda environment for DESeq2  data analysis 
with R and bioconductor in my window system using Windows
Subsystem Linux (WSL)

## Code in WSL terminal
- start the WSL in windows
- Type the following to create a conda environment for R4.5 name as R4.5WS
```{cmd}
conda create --name R4.5WS
```
- Activate the R4.5WS conda environment as follows
```{CMD}
conda activate R4.5WS
```
- Install R4.5 in the R4.5WS conda environment
```{CMD}
conda install conda-forge::r-base
```
`conda-forge` is the channel from anaconda.org. [conda-forge channel for r-base 4.5.2](https://anaconda.org/conda-forge/r-base)
