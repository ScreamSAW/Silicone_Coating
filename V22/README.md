# V22 LAMMPS data generator

`v22_generator.cpp` creates a LAMMPS `atom_style full` data file containing coordinates, bonds, angles, and dihedrals for the four-component coarse-grained V22 formulation.

1. `N1`, `M1`: bifunctional network strands. The two terminal beads are reactive atom type 2.
2. `N2`, `M2`: multifunctional cross-linkers. Reactive sites are atom type 3.
3. `N3`, `M3`: neutral linear PDMS oil filler. All beads are atom type 1.
4. `N4`, `M4`: five-bead star-like moderators. The center is neutral type 1 and the four arms are reactive type 2.

Atom type 1 is neutral. Types 2 and 3 may react during later simulations. The source preserves the legacy declaration of four atom types, although atom type 4 is currently unused.

## Build

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic v22_generator.cpp -o v22_generator
```

Display the built-in command reference with:

```bash
./v22_generator --help
```

## Default setting

Run the default formulation with:

```bash
./v22_generator
```

| Component | Meaning | Chain size `N` | Molecules `M` | Bead types |
|---|---|---:|---:|---|
| 1 | Bifunctional network strand | `N1=128` | `M1=900` | Type 1 interior; type 2 ends |
| 2 | Cross-linker | `N2=32` | `M2=225` | Type 1 neutral; type 3 reactive sites |
| 3 | Linear PDMS filler | `N3=32` | `M3=425` | Type 1 |
| 4 | Star-like moderator | `N4=5` | `M4=6` | Type 1 center; type 2 arms |

Other defaults:

| Setting | Default |
|---|---:|
| Cross-linker functionality | `8` |
| Cross-linker site distribution | `random` |
| Cross-linker random seed | `20260722` |
| PDMS loading | approximately `10 wt%` |
| Bead mass | `74 g/mol` |
| Target density | `0.1 g/cm^3` |
| Initial bond length | `2.801 angstrom` |
| Molecule placement spacing | `7.0 angstrom` |
| Geometry | Cubic bulk box |
| Star-moderator coordinate seed | `5489` |
| Output filename | `data.V22_PDMS_N32_10wt` |

`M2=225` is derived from stoichiometry, not independently specified. `M3=425` gives approximately 10 wt% PDMS after integer rounding.

## All command-line inputs

| Option | Value | Default | Meaning and restrictions |
|---|---|---:|---|
| `--n1` | integer | `128` | Beads per network strand; nonnegative |
| `--m1` | integer | `900` | Number of network strands; determines stoichiometric `M2` |
| `--n2` | integer | `32` | Beads per cross-linker; must be at least the functionality |
| `--n3` | integer | `32` | Manual PDMS chain length |
| `--m3` | integer | `425` | Manual PDMS chain count; overridden by `--filler-wt` |
| `--n4` | integer | `5` | Moderator size; must remain 5 when `M4>0` |
| `--m4` | integer | `6` | Number of star-like moderators |
| `--filler-length` | integer | `32` | Recommended PDMS chain-length input; alias for `--n3` |
| `--filler-wt` | number | `10` conceptually | PDMS percent of total weight; range `0 <= X < 100`; calculates `M3` |
| `--functionality` | integer | `8` | Reactive sites per cross-linker; `3 <= F <= min(16,N2)` |
| `--crosslink-distribution` | text | `random` | Either `random` or `regular` |
| `--crosslink-seed` | integer | `20260722` | Reproducible random type-3 site selection |
| `--mass` | number | `74` | Common bead mass in g/mol; must be positive |
| `--density` | number | `0.1` | Target density in g/cm3; must be positive |
| `--bond-length` | number | `2.801` | Initial bond length in angstrom; must be positive |
| `--spacing` | number | `7.0` | Initial molecule-placement spacing; must be positive |
| `--thickness` | number | omitted | Enables film mode with fixed `Lz` in angstrom |
| `--seed` | integer | `5489` | Reproducible star-moderator coordinate seed |
| `--output` | path | automatic | Overrides automatic output naming |
| `--help` | none | — | Prints built-in help and exits |

`--m2` is deliberately rejected because the program calculates `M2` from stoichiometry.

For normal V22 work, leave components 1, 2, and 4 fixed and vary the PDMS using `--filler-length` and `--filler-wt`. The lower-level `--n3` and `--m3` inputs are retained for manual studies.

## PDMS filler controls

Change the PDMS chain length and weight percentage with:

```bash
./v22_generator --filler-length 64 --filler-wt 15
```

The V22 base consists of components 1, 2, and 4. With equal bead masses, the program calculates:

```text
base_beads   = N1*M1 + N2*M2 + N4*M4
filler_beads = base_beads * (filler_wt/100) / (1 - filler_wt/100)
M3           = round(filler_beads / filler_length)
```

The example writes `data.V22_PDMS_N64_15wt`. Since `M3` must be an integer, the program prints the realized weight percentage, which may differ slightly from the request.

Either input may be used alone:

```bash
./v22_generator --filler-length 64  # retains 10 wt%
./v22_generator --filler-wt 15      # retains N32
./v22_generator --filler-wt 0       # no PDMS filler
```

Use `--output FILE` to replace the automatically generated filename.

## Film geometry

Omitting `--thickness` preserves the original cubic bulk system. Supplying a
positive thickness enables film mode:

```bash
./v22_generator --thickness 100
```

The requested value becomes the fixed `Lz`. The program preserves the target
density by calculating equal lateral dimensions:

```text
volume = total_bead_mass / density
Lx = Ly = sqrt(volume / Lz)
```

The example automatically writes
`data.V22_PDMS_N32_10wt_film_Lz100`. Film geometry can be combined with the
other composition controls:

```bash
./v22_generator \
  --thickness 80 \
  --filler-length 64 \
  --filler-wt 15
```

The LAMMPS data file defines the rectangular simulation box but does not define
boundary behavior. To keep the film thickness fixed during simulation, use
nonperiodic boundaries in `z` and avoid barostatting `z`. For example, the
LAMMPS input script may use `boundary p p f` together with an appropriate
`fix wall/*` command at the lower and upper z boundaries. The exact wall style
and parameters depend on the desired physical surface interaction.

## Cross-linker stoichiometry

A network requires cross-linker functionality of at least 3. The maximum supported functionality is 16, and functionality cannot exceed `N2`.

Every network strand has two reactive type-2 ends, so the program enforces:

```text
M2 * functionality = 2 * M1
M2 = 2*M1/functionality
```

Exact stoichiometry requires integer `M2`. If `2*M1` is not divisible by the requested functionality, generation stops with an error rather than silently creating an imbalance.

For the default `M1=900`, the exactly compatible functionalities between 3 and 16 are:

```text
3, 4, 5, 6, 8, 9, 10, 12, 15
```

Examples:

```bash
./v22_generator --functionality 8   # M2 = 225
./v22_generator --functionality 12  # M2 = 150
./v22_generator --functionality 16  # error: 1800/16 is not an integer
```

## Reactive-site distribution

Regular distribution places type-3 sites at beads `1,3,5,7,...`:

```bash
./v22_generator \
  --functionality 8 \
  --crosslink-distribution regular
```

Regular placement additionally requires `functionality <= ceil(N2/2)`.

Random distribution selects distinct positions from `1...N2`:

```bash
./v22_generator \
  --functionality 12 \
  --crosslink-distribution random \
  --crosslink-seed 314159
```

The same cross-linker seed reproduces the same positions. `--crosslink-seed` affects only reactive-site selection; `--seed` independently controls star-moderator coordinates. Selected cross-linker sites are printed during generation.

## Complete example

This example retains the fixed V22 network strands and moderators, uses functionality-12 cross-linkers with random sites, and adds N64 PDMS at 15 wt%:

```bash
./v22_generator \
  --functionality 12 \
  --crosslink-distribution random \
  --crosslink-seed 314159 \
  --filler-length 64 \
  --filler-wt 15 \
  --seed 5489
```

The program calculates `M2=150`, calculates the required integer `M3`, reports the realized composition, and writes `data.V22_PDMS_N64_15wt`.

## Modeling notes

- Topology counts in the LAMMPS header are derived from the topology actually generated.
- Linear molecules use bond type 1. Star-moderator arms use bond type 2.
- The star moderator preserves the original planar five-bead geometry and its six arm-arm angles.
- Star-moderator centers preserve the original narrow z slab, from `0.02L` to `0.12L`.
- Coordinates are initial placements and can overlap under periodic wrapping. Energy minimization or another packing/overlap-removal procedure is needed before production dynamics.
- The filler calculation currently assumes the same coarse-grained bead mass for all components.
