AMELLIE Macros
===================

Contributors: Charlie, Lisa

## Usage

#### runSimulation.mac

An basic example macro to simulate an AMELLIE run.

Requires a few changes to run:
1. A valid GEO file must be used
2. A valid LED must be used
3. A valid FIBRE must be used
4. An output file must be specified


#### convertZDAB.mac

This macro will take a raw ZDAB file and apply calibration, giving a ROOT file as the output. Simply run this macro with the ```-i``` and ```-o``` flags, e.g:

```rat convertZDAB.mac -i /path/to/your/file.zdab -o /path/to/your/file.root```