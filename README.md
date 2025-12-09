<a name="top"></a>

# Overset-Exploded

> Overse-Exploded, overset post-processor, automatic blocks-splitting, load-balancing, blocks-explosion

| [Main Features](#main-features) | [Usage](#usage) | [Compile](#compile) | [Download](#download) | [Authors](#authors) | [Copyrights](#copyrights) |

## Main Features

Overset-Exploded is an overset post-processor: it takes in input overset outputs and generate a new format
per-block inputs for Xnavis/Xall. In particular, the overset grd and icc outputs are loaded and automatically
split in order to load-balance the global workload accordingly to the number of processes used. The global rcc
array is split per-block and the new adjacent data for split blocks are added. It also generate a proc.input
file with the load-balanced distribution. The level of admissible unbalancing can be customized in input, as
well as the multigrid level to be preserved when splitting blocks.

Go to [Top](#top)

## Usage

Overset-Exploded has a built in help manual:

```bash
overset-exploded: overset post-processor, automatic blocks-splitting, load-balancing, blocks-explosion
usage:
   overset-exploded [args]
args list:
   -grd file_name_grd               => GRD file name, default "cc.01.grd"
   -icc file_name_icc               => ICC file name, default "cc.01"
   -proc-input file_name_proc_input => proc.input file name, default "proc.input"
   -np processes_number             => number of processes for load balancing, default 1
   -max-unbalance mu                => maximum processes unbalancing in percent, default 1%
   -mgl mgl                         => multigrid level to be preserved, default 4
   -tec                             => enable tecplot output for debug, default .false.
   -save-imploded                   => save imploded blocks after explosion, default .false.
   -save-exploded                   => save exploded blocks, default .false.
   -exploded-basename               => exploded files basename, default "exploded-"
   -h, --help                       => print this help message
examples:
   overset-exploded -np 32
   overset-exploded -grd cc.02.grd -icc cc.02 -np 16
   overset-exploded -np 16 -max-unbalance 4
   overset-exploded -np 16 -proc-input proc.input-pes16
   overset-exploded -np 16 -proc-input proc.input-pes16 -save-imploded
   overset-exploded -np 16 -proc-input proc.input-pes16 -save-exploded
   overset-exploded -grd cc.03.grd -icc cc.03 -np 2 -tec -max-unbalance 3

```

Typical usage:

+ `overset-exploded -np 32` balance the mesh for 32 processes, split blocks (if possible) until maximum unbalancing is below 1%;
  only proc.input is saved, split blocks are not saved, neither in old imploded-legacy format nor in new exploded one;
+ `overset-exploded -np 16 -max-unbalance 4` the same as before, but the maximum unbalancing is admitted is 4%;
+ `overset-exploded -np 16 -proc-input proc.input-pes16 -save-imploded` the load-balancing is done for 16 processes, the
  proc.input is save as `proc.input-pes16` and the split-balanced blocks are save in old imploded format;

Go to [Top](#top)

## Compile

Overset-Exploded is a self-contained, single file, Fortran code that can be easily compiled with any Fortran compiler supporting
ISO 2008 standard. A makefile and fobos files are provided to build it by means of GNU make and FoBiS.py builder.

### Compile with GNU Make

In the root of the repository type

```bash
make
```

The program will be assembled in `exe/overset-exploded`.

### Compile with FoBiS.py

In the root of the repository type

```bash
FoBiS.py build
```

The program will be assembled in `exe/overset-exploded`.

Go to [Top](#top)

## Download

Use git to clone this repository:

```bash
git clone https://github.com/szaghi/overset-exploded
```

Go to [Top](#top)

## Authors

+ Andrea di Mascio, [andrea.dimascio@univaq.it](andrea.dimascio@univaq.it)
+ Stefano Zaghi, [stefano.zaghi@cnr.it](stefano.zaghi@cnr.it)

Go to [Top](#top)

## Copyrights

Released under GPL v3, a copy of which is included into this repository.

Go to [Top](#top)
