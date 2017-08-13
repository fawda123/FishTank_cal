## README

Contents include FishTank executable for Linux and R functions for calibration.

### Compile FishTank in CGEM repo

From `~/CGEM` repository, compile the latest version of FishTank:
```bash
cd CGEM
module add intel
source modules.sh
make
```

The new compiled version of FishTank can be copied to a new directory:

```bash
cp CGEM/FishTank.exe ../FishTank_cal
```

Available modules can be viewed:
```bash
module avail
```

Added modules can be viewed:
```bash
module list
```

The CGEM_InputFile lives in the CGEM folder.  All FishTank files are in `CGEM/data/FishTank`.  The file `CGEM/data/MyFiles.inp` must have the full path to the FishTank folder.

FishTank is executed:
```bash
./FishTank.exe
```

The output file is in `~/CGEM/NETCDF`.  It can be renamed and R plots can be created:
```bash
cd NETCDF
mv output.000000.nc cgem_0D.nc
Rscript make_plots_0D.R cgem
evince cgem_0D.pdf
```

### Running FishTank in current repo

The appropriate modules must be loaded:

```bash
module add intel
source ../CGEM/CGEM/modules.sh
./FishTank.exe
```
