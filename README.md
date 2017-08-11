## README

Contents include FishTank executable for Linux and R functions for calibration.

FishTank executable was created from `~/CGEM` repository at http://github.com/USEPA/CGEM, master branch.  Always track development in the CGEM branch.

From `~/CGEM` repository, compile the latest version of FishTank:
```bash
cd CGEM
module add intel
source modules.sh
make
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
