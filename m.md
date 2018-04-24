**sDMD**

**User Manual**

![](media/image1.png){width="3.5075240594925634in" height="3.5in"}

*Version 0.93*

Table of Contents {#table-of-contents .TOCHeading}
=================

[1. Introduction 1](#introduction)

[2. Units 4](#units)

[3. Models and Interactions 5](#models-and-interactions)

[3.1 Models 5](#models)

[3.1.1 Intermediate-resolution Coarse Grained Model
5](#intermediate-resolution-coarse-grained-model)

[3.1.2 High-resolution All-atom Model
10](#high-resolution-all-atom-model)

[3.2 Interactions 13](#interactions)

[3.2.1 Algorithms 13](#algorithms)

[3.2.2 Energy Minimization 17](#energy-minimization)

[3.2.3 Confined Walls 18](#confined-walls)

[3.2.4 Replica Exchange Molecular Dynamics Simulation (REMD)
20](#replica-exchange-molecular-dynamics-simulation-remd)

[4. Program Structure 23](#program-structure)

[5. Program Usage 25](#program-usage)

[5.1 Installation 25](#installation)

[5.2 Run Parameters 25](#run-parameters)

[5.2.1 Non-REMD simulation 26](#non-remd-simulation)

[5.2.2 REMD Simulation 28](#remd-simulation)

[5.2.3 Analysis 30](#analysis)

[6. Format 32](#format)

[6.1 Coordinate File 32](#coordinate-file)

[6.1.1 GRO File 32](#gro-file)

[6.1.2 PDB File 34](#pdb-file)

[6.2 Force Field File 35](#force-field-file)

[6.2.1 Parameters of amino acid 35](#parameters-of-amino-acid)

[6.2.2 Potential Table 37](#potential-table)

[7. Functions 38](#functions)

[7.1 Simulation Source Code 38](#simulation-source-code)

[DMD.c 38](#dmd.c)

[DataSave.c 38](#datasave.c)

[Event.c 38](#event.c)

[Initialization.c 39](#initialization.c)

[List.c 39](#list.c)

[Models.c 39](#models.c)

[REMD.c 39](#remd.c)

[SGThread.c 39](#sgthread.c)

[ThreadProcess.c 39](#threadprocess.c)

[TimePrediction.c 40](#timeprediction.c)

[ToolFunctions.c 40](#toolfunctions.c)

[sREMD.c 40](#sremd.c)

[sServer.c 40](#sserver.c)

[7.2 Analysis Source Code 40](#analysis-source-code)

[Analysis.c 40](#analysis.c)

[FileManage.c 41](#filemanage.c)

[SystemInformation.c 41](#systeminformation.c)

[Cluster.c 41](#cluster.c)

[Energy.c 41](#energy.c)

[HBRamach.c 41](#hbramach.c)

[PBCAdjust.c 41](#pbcadjust.c)

[REMD.c 41](#remd.c-1)

[Tools.c 41](#tools.c)

**Introduction**
================

Discrete molecular dynamics (DMD) simulation algorithm is an extremely
fast alternative to traditional molecular dynamics. It was first
introduced in 1959 by Alder and Wainwright^1^ for simulations of hard
spheres. Later it was used by Rapaport^2-4^ for simulation of polymer
chains. Nowadays it was adopted for simulations of protein-like
polymers. It is extremely fast and suitable for the simulations of large
systems on long time scales.

*Figure 1-1 DMD potential (a) hard sphere collision (b) attractive
square well interaction (c) repulsive soft sphere interaction (d)
covalent bond and auxiliary bond interaction*

During DMD simulation, potentials applied on particles are approximated
by discontinuous step-functions of inter-particle distance $r$ (see
Figure *1-1*). Thus no forces will be exerted on particles until their
distance becomes equal to the point of a discontinuity on the potential.
When particle encounters a potential discontinuity, we call it have an
"event". Between events, particles move at constant velocities. Because
the energy change during each event is known, the post-event velocities
can be calculated by solving the conservation of momentum and the
conservation of energy simultaneously. Thus, the trajectory of a
particle can be simulated discontinuously between events. The simulation
code for DMD is significantly different than that for traditional
molecular dynamics because the particles are moved discontinuously from
one event to the next with known velocities.

Here we proudly introduce the sDMD, a simulation package based on the
DMD technique and a high-resolution all-atom molecular model. The
package is written by C language from scratch and has been highly
optimized. Its modularized structure has it be highly expandable: new
molecular groups can be introduced by simply adding new formatted
connection and potential data files, and no modification on the code
would be necessary. The package, including the simulation executable and
analysis executable, is also very easy to use: the configurations are
controlled by an easy-to-read configuration file and several command
flags.

Currently the package is capable to run simulations under various
situations as follows,

-   NVT and NVE ensembles

-   bulk conditions as well as parallel pores, cylindrical tubes and
    spherical cavities

-   replica exchange molecular dynamics simulation (REMD)

-   dynamic confined systems as the confined walls expand, compress and
    pulse.

-   molecule flow through a tube between two reservoirs

and the analysis executable can be used to analyze the following,

-   remove periodic boundary conditions of trajectory

-   compute the number of different types of hydrogen bonds

-   compute Ramachandran plot

-   compute system temperature, kinetic energy and potential energy

-   compute number of clusters during aggregation

-   analyze data from T-REMD by weighted histogram analysis method
    (WHAM)

This manual will describe the theories behind the DMD simulation, as
well as the models of molecules and algorithms of interactions applied
in the MD simulation. It will also describe the function of each file in
the source code and demonstrate the ways to install and use the package.

1\. Alder, B. J.; Wainwright, T. E., Studies in Molecular Dynamics. I.
General Method. *The Journal of Chemical Physics* **1959,** *31* (2),
459-466.2. Rapaport, D. C., Molecular dynamics simulation of polymer
chains with excluded volume. *Journal of Physics A: Mathematical and
General* **1978,** *11* (8), L213.3. Rapaport, D. C., Molecular dynamics
study of a polymer chain in solution. *The Journal of Chemical Physics*
**1979,** *71* (8), 3299-3303.4. Rapaport, D., The event scheduling
problem in molecular dynamic simulation. *Journal of Computational
Physics* **1980,** *34* (2), 184-201.

**Units**
=========

*Table 2-1 Basic units*

+------------+----------+---------------------------------+
| *Quantity* | *Symbol* | *Unit*                          |
+============+==========+=================================+
| length     | L        | Å = 10^-10^ m                   |
+------------+----------+---------------------------------+
| mass       | m        | u (atomic mass unit)            |
|            |          |                                 |
|            |          | = 1.6605402(10) x 10^-27^ kg    |
|            |          |                                 |
|            |          | (1/12 the mass of a ^12^C atom) |
+------------+----------+---------------------------------+
| energy     | E        | kcal mol^-1^                    |
|            |          |                                 |
|            |          | = 4.184 kJ mol^-1^              |
|            |          |                                 |
|            |          | = 4184 m^2^ kg s^-1^ mol^-1^    |
+------------+----------+---------------------------------+

*Table 2-2 Derived units*

  *Quantity*    *Symbol*   *Unit*
  ------------- ---------- -----------------------------------
  time          T          50 fs = 50 x 10^-15^ s
  velocity      V          Å (50 fs)^-1^ = 2 x 10^3^ m s^-1^
  temperature   T          503 K
  force         F          u Å (50 fs)^-2^ = 66.4 pN

*Table 2-3 Constants*

  *Symbol*   *Name*                 *Unit*
  ---------- ---------------------- ---------------------------------------
  *N~A~*     Avogadro's number      6.0221367(36) x 10 23 mol^-1^
  *k~B~*     Boltzmann's constant   1.987204118 x 10-3 kcal mol^-1^ K^-1^

**Models and Interactions**
===========================

This section will describe the molecular models and interaction
algorithms that used in the sDMD.

**Models**
----------

Different protein models can be used in DMD ranging from high-resolution
all-atom model to low-resolution coarse-grained model. The latest
version of the sDMD only uses the high-resolution all-atom model. The
early version uses coarse grained model. Both models share a lot of
characteristics. Here will give a description of the models starting
from the low-resolution and then upgrade to the high-resolution.

### **Intermediate-resolution Coarse Grained Model**

In the early version, the sDMD uses an intermediate-resolution coarse
grained model, named PRIME, developed by Smith and Hall^1-4^. It
correctly reproduces the backbone geometry and has been highly
successful in reproducing several important and experimentally
determined features of the proteins under bulk conditions. In RPIME,
every amino acid is represented by four beads: N, C~α~, C, and R (see
*Figure 3-1*). The N bead represents the amide nitrogen and its
hydrogen; the C~α~ bead represents the alpha-carbon and its hydrogen;
the C bead represents the carbonyl carbon and oxygen; the R bead
represent the side chain, all of which are assumed to have the same
diameter as CH~3~ (alanine). The inter-peptide bond is assumed to be in
the *trans* configurations, all the backbone bonds' lengths and bond
angles are fixed at their experimentally measure values, and the
distance between consecutive C~α~ bead is also fixed according to
experimental data. *Table 3-1* presents all the relevant parameters of
the model. Solvent molecules are not explicitly included in the model.
The effect of solvent is factored into the energy function as a
potential of mean force.

![](media/image2.png){width="5.021738845144357in" height="2.0in"}

*Figure 3-1 The four-bead model of the protein backbone. Bold lines show
covalent bond, dashed lies show auxiliary bonds which helps maintain
correct backbone geometry, and thick broken lines show possible hydrogen
bonds.^5^*

*Table 3-1 Parameter values for the proteins used in the DMD
simulations*

  **Beads**               **Diameter (Å)**
  ----------------------- ------------------
  N                       3.300
  C~α~                    3.700
  C                       4.000
  R                       4.408
                          
  **Bonds**               **Length (Å)**
  N~i~-C~α,i~             1.460
  C~α,i~-N~i+1~           1.510
  C~i~-N~i+1~             1.330
  C~α,i~-R~i~             1.531
                          
  **Auxiliary Bonds**     **Length (Å)**
  N~i~-C~i~               2.450
  C~α,i~-N~i+1~           2.410
  C~i~-C~α,i+1~           2.450
  C~α,i~-C~a,i+1~         3.800
  N~i~-R~i~               2.440
  C~i~-R~i~               2.490
                          
  **Bond Angles**         **Angle (deg)**
  ∠N~i~-C~α,i~-C~i~       111.0
  ∠C~α,i~-C~i~-N~i+1~     116.0
  ∠C~i~-N~i+1~-C~α,i+1~   122.0
  ∠R~i~-C~α,i~-C~i~       109.6
  ∠R~i~-C~α,i~-C~i~       110.1

Beads in the protein are subject to five different types of forces
during events: (1) infinite repulsion due to excluded volume effect
(hard-sphere collision), (2) infinite attraction between the beads
connected with covalent bonds and auxiliary bonds, (3) attraction
between pairs of back-bone beads during hydrogen bond formation, (4)
attraction between hydrophobic side chains, and (5) repulsion between
the neighborhood beads of hydrogen bond (soft-sphere collision). All
these forces can be represented by the following step-functions,

  -- ---------------------------------------------------------- -------------
     $$U_{\text{ij}}\left( r \right) = \left\{ \begin{matrix}   *(Eqn 3-1)*
     \infty\ \ \ \ \ \ \ \ \ \ r \leq \sigma_{\text{ij}} \\     
     0\ \ \ \ \ \ \ \ \ \ \ r > \sigma_{\text{ij}} \\           
     \end{matrix} \right.\ $$                                   
  -- ---------------------------------------------------------- -------------

This is for a hard-sphere potential. $r$ is the distance between beads
$i$ and $j$; $\sigma_{\text{ij}}$ is the sum of bead radius.

  -- ---------------------------------------------------------------------------------------- -------------
     $$U_{\text{ij}}\left( r \right) = \left\{ \begin{matrix}                                 *(Eqn 3-2)*
     \infty\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ r \leq \sigma_{\text{ij}} \\     
      - \varepsilon\ \ \ \ \ \ \ \ \ \ \sigma_{\text{ij}} < r \leq \text{λσ}_{\text{ij}} \\   
     0\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ r > \text{λσ}_{\text{ij}} \\            
     \end{matrix} \right.\ $$                                                                 
  -- ---------------------------------------------------------------------------------------- -------------

This is for a square-well or shoulder potential depending on the sign of
$\varepsilon$ (positive for well depth, attractive; negative for
shoulder height, repulsive). $\text{λσ}_{\text{ij}}$ is the
well/shoulder diameter.

  -- -------------------------------------------------------------------------------------------------------- -------------
     $$U_{\text{ij}}\left( r \right) = \left\{ \begin{matrix}                                                 *(Eqn 3-3)*
     \infty\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ r \leq l\left( 1 - \delta \right) \\   
     0\ \ \ \ \ \ \ \ \ \ l\left( 1 - \delta \right) < r < l\left( 1 + \delta \right) \\                      
     \infty\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ r \geq l\left( 1 + \delta \right) \\   
     \end{matrix} \right.\ $$                                                                                 
  -- -------------------------------------------------------------------------------------------------------- -------------

This is for a bond interaction. $l$ is the bond length and $\delta$ is
the bond vibration tolerance.

Note that, as explained in the work by Takada et al,^6^ interaction
involving a pair of beads that are very close on the chain are more
faithfully represented by an interaction between the atoms themselves,
not the united atoms or the coarse grained beads. In other words, for
example, neighboring carbonyl carbon atoms experience interactions based
on the diameter of a carbon atom, not on the diameter of the C bead.
Consequently, for interaction between pairs of beads that separated by
three or fewer bonds, the beads are allowed to overlap in their
simulations by up to 25% of their diameters; for those separated by four
bonds, they are allowed to overlap up to 15%.

Backbone beads and C~α~ and R beads are connected by covalent bonds.
Bonds are allowed to move freely over a small range between
$l\left( 1 - \delta \right)$ and $l\left( 1 - \delta \right)$. The
choice of $\delta$ defines the acceptable range of bond vibration. In
the simulation, $\delta$ is chosen to be 0.02375. Auxiliary bonds are
used between next-neighbor beads along the backbone of the chain to hold
backbone angles fixed. For example, an auxiliary is placed between N~i~
and C~i~ in *Figure (3-1)* to force N~i~-C~a,i~-C~i~ angle to be near
its ideal value. Auxiliary bonds are also included between neighboring
pairs of C~α~ beads to maintain their distances near the experimentally
determined constant value. The combination of C~a,i~-N~i+1~,
C~i~-C~a,i+1~, and C~a,i~-C~a,i+1~ auxiliary bonds and the C~a,i~-C~i~,
C~i~-N~i+1~, and N~i+1~-C~a,i+1~ covalent bonds holds the inter-peptide
group (C~a,i~-C~i~-N~i+1~-C~a,i+1~) in the *trans* configuration.
Auxiliary bonds are also placed between side-chains and N and C beads on
backbone. This is to hold the side-chain beads fixed relative to the
backbone so that all model residues are *L*-isomers.

Hydrogen bond (HB) can be formed between N and C beads on backbones.
Four criteria have to be met (see *Figure 3-2*): (1) the distance
between N~i~ and C~j~ beads is 4.2 Å; (2) the neighbor beads, C~i-1~ and
C~a,i~, C~a,j~ and N~j+1~, are at appropriate positions; (3) neither N
nor C bead is already involved in a hydrogen bond; (4) N and C beads are
separated by at least three intervening residues or they are in
different peptides. To ensure that criteria (2) is satisfied, we require
that the four atom pairs N~i~-C~a,j~, N~i~-N~j+1~, C~j~-C~a,i~,
C~j~-C~i-1~ shown in *Figure (3-2)* be separated by a distance greater
than $d_{\text{ij}}$, which is chosen to maintain the hydrogen bond
angle constraints. Their values are given in *Table 3-2*.

When two beads N~i~ and C~j~ approach each other at the distance
$d_{\text{ij}} = 4.2\mathring{\mathrm{A}}$, we will evaluate the total
potential energy change by checking the four auxiliary interactions. The
potential energy change can be $- \varepsilon_{\mathrm{\text{HB}}}$,
$0$, $\varepsilon_{\mathrm{\text{HB}}}$,
$2\varepsilon_{\mathrm{\text{HB}}}$ and
$3\varepsilon_{\mathrm{\text{HB}}}$, depending on the orientation of the
N~i~, C~j~, and their neighbors. If the kinetic energy is enough to
overcome the total potential energy change, the forward reaction will
happen. N~i~ and C~j~ beads will change their type into N~i~' and C~j~',
respectively, and interact with each other (the other interactions with
N~i~ and C~j~ will not be affected) according to the interaction
parameters related to their new types. Otherwise, the two atoms N~i~ and
C~j~ do not change their types and undergo original hard-sphere
collision. The thermal fluctuation distort the orientation of the
hydrogen bond and large fluctuations may break the hydrogen bond. Once
the two atoms N~i~' and C~j~' come again to the exact distance
$d_{\text{ij}} = 4.2\mathring{\mathrm{A}}$, a reverse reaction may
happen. We will check the potential energy change by using the similar
way we check in HB formation, and decide if the HB will break so that
N~i~' and C~j~' will change back to their original types, N~i~ and C~j~,
respectively.

![C:\\Users\\Size\\Desktop\\Untitled-1.jpg](media/image3.jpeg){width="2.0616765091863516in"
height="3.0in"}

*Figure 3-2 Hydrogen bond structure of coarse grained model. Thick
dashed lines indicate hydrogen bond; thin dashed lines indicate
auxiliary bonds helping to maintain hydrogen bond orientation, and bold
lines show the covalent bonds of the backbone.*

*Table 3-2 Parameter values for hydrogen bond*

  **Pairs**   $$\mathbf{d}_{\mathbf{\text{ij}}}\left( \mathbf{\mathring{\mathrm{A}}} \right)$$
  ----------- ----------------------------------------------------------------------------------
  Ni-Cα,j     5.00
  Ni-Nj+1     4.74
  Cj-Cα,i     4.86
  Cj-Ci-1     4.83

To maintain the temperature around the target value in NVT ensemble, we
implement Andersen thermostat.^7^ With this procedure, all beads in the
simulation are subjected to random, infrequent collisions with ghost
particles having the same mass as themselves. The post-event velocity of
a bead colliding with a ghost particle is chosen randomly from a
Maxwell-Boltzmann distribution centered at the simulation temperature.
The collision frequency is controlled at 1% of the total events number.

The solvent molecules are not explicitly included. They are modeled
implicitly by means of a square well attraction (potential of mean
force) between two hydrophobic side-chain beads. The side-chain beads, R
beads, in a peptide are assigned as either "H" or "P": "H" means it is
hydrophobic and such R beads will attract with each other when their
distance is equal to $1.5\sigma_{\text{ij}}$ and as long as they are
separated by at least three intervening residues in the same peptide or
they are on the different peptides; "P" means the side-chain bead is
polar. The interactions between two polar beads or between a polar and a
hydrophobic side-chain beads will be a normal hard sphere collision.

### **High-resolution All-atom Model**

In order to apply the sDMD to more complex systems, we upgraded the
molecular and potential models from coarse grained to all-atom. Large
amounts of codes were rewritten and optimized to fulfill the more
modularized and efficient architecture. Details of the format of the
library files will be described in Section 5.

The all-atom models are developed by Ding et al^8^, and thanks for the
helps of Dr. Ding, we got the updated force field data (the data on the
paper was outdated). The basic concepts and algorithms used in the
all-atom models are similar to those in the coarse grained models: it
still uses auxiliary bonds to maintain the configuration of molecules,
and it still uses the step functions to represent the forces. But now,
after the updates, the model represents all the atoms explicitly (except
for some light hydrogens), see Figure 3-3; the step function is more
realistic and closer to the continuous force, see Figure 3-4(a). Due to
the updated multi-step potential profile, now most interactions between
atoms combine both attraction and repulsion, unlike the previous that
only has either one of them. It also adds a new type of force, the
dihedral constraint. Although it is maintained by an auxiliary bond, it
has more complex profile than a regular bond, see Figure 3-4(b).

![](media/image4.JPG){width="3.48954615048119in" height="3.0in"}

*Figure 3-3 Schematic diagram for the all-atom protein model. The solid
lines represent the covalent bonds and the dashed lines represent the
auxiliary bonds.*

![](media/image5.JPG){width="3.5in"
height="2.4270548993875765in"}![](media/image6.png){width="2.992788713910761in"
height="2.5in"}

*Figure 3-4 (a) Non-bonded interactions in all-atom DMD. The continuous
dashed line corresponds to the VDW and solvation interaction between two
carbon atoms. The solid step function is the discretized potential for
DMD. (b) Potential profile of dihedral constraint.*

Since the all-atom model represents the oxygen and hydrogen atoms
explicitly, the structure of hydrogen bond is also optimized. The new
hydrogen bond structure is shown in Figure 3-5.

![C:\\Users\\Size\\AppData\\Local\\Microsoft\\Windows\\INetCache\\Content.Word\\HB.JPG](media/image7.jpeg){width="4.5in"
height="1.752584208223972in"}

*Figure 3-5 Hydrogen bond structure of all-atom model. Hi is the
hydrogen atom; Aj is the acceptor and could be different types of atom;
Di is the donor and Xj is the heavy atom directly bonded to Aj. The
thick solid lines indicate the covalent bonds between the atoms; the
thin solid lines indicate the hydrogen bond; the dashed lines indicate
the auxiliary bonds for helping maintain the hydrogen bond orientation.*

By adjusting the type of Aj atom, the hydrogen bond can be formed
between backbone-backbone, backbone-side chain, and side chain-side
chain. The donors include backbone amide hydrogen atoms and side chain
polar hydrogen atoms of His, Trp, Tyr, Asn, Gln, Arg, and Lys. The
acceptors include backbone carbonyl oxygens, side chain oxygen of Asp,
Glu, Ser, Thr, and Tyr, and the side chain nitrogen of His.

The mechanism to model the hydrogen bond forming and breaking is similar
to the coarse grained model: once the hydrogen, Hi, and the acceptor,
Aj, reach the interaction, the distances between HiXj and DiAj, which
define the orientation of the hydrogen bond, will be evaluated; the
total potential energy change, ΔE, between Hi/Aj and other surrounding
atoms are also evaluated before and after the putative hydrogen bond
formation; if these distances satisfy the preset range, and the total
kinetic energy is enough to overcome the potential energy change, ΔE,
the hydrogen bond will be allowed to form, and forbad otherwise; after
the formation of hydrogen bond, the atoms Hi and Aj will become Hi' and
Aj', respectively, and they will interact with other atoms by using the
new pairs of potential.

**Interactions**
----------------

Here will describe the basic algorithms behind the DMD simulation. It
will also give a talk about the mechanisms in energy minimization, wall
interaction, and REMD.

### **Algorithms**

This appendix will describe how DMD works.

Assume we have two particles, $i$ and $j$, each with diameter σ. The
coordinates of the particles are the vector
${\overrightarrow{r}}_{i} = \left\lbrack x_{i},y_{i},z_{i} \right\rbrack$
and
${\overrightarrow{r}}_{j} = \left\lbrack x_{j},y_{j},z_{j} \right\rbrack$
and they have velocity vector ${\overrightarrow{v}}_{i}$ and
${\overrightarrow{v}}_{j}$. The collision event and collision time can
be determined by using the relative position and relative velocity,
${\overrightarrow{r}}_{\text{ij}}$ and
${\overrightarrow{v}}_{\text{ij}}$.

*Figure 3-6 Sample System (2D)*

A collision event is dependent of whether the distance between particles
$i$ and $j$ will become equal to the average molecular dimeter
$\sigma_{\text{ij}} = \left( \sigma_{i} + \sigma_{j} \right)/2$ sometime
in the future. Using superscript "o" to represent an initial condition,
the relative position vector at any time in the future is found by the
following equation,

  -- -------------------------------------------------------------------------------------------------------------------------------------------- -----------
     $${\overrightarrow{r}}_{\text{ij}} = {\overrightarrow{r}}_{\text{ij}}^{o} + {\overrightarrow{v}}_{\text{ij}}^{o}\left( t - t^{o} \right)$$   (Eqn 3-4)
  -- -------------------------------------------------------------------------------------------------------------------------------------------- -----------

So if ${\overrightarrow{r}}_{\text{ij}}$ and
${\overrightarrow{v}}_{\text{ij}}$ have opposite sense, the particles
will be approaching and a collision is possible. We may check to see if
the dot product
${\overrightarrow{r}}_{\text{ij}} \bullet {\overrightarrow{r}}_{\text{ij}}$
becomes equal to $\sigma_{\text{ij}}^{2}$ in the future at some time
increment $\left( t - t^{o} \right)$,

     $$\sigma_{\text{ij}}^{2} = {\overrightarrow{r}}_{\text{ij}} \bullet {\overrightarrow{r}}_{\text{ij}}$$                                                                                                                                                                                                                      (Eqn 3-5)
  -- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- -----------
     $$\sigma_{\text{ij}}^{2} = \left\lbrack {\overrightarrow{r}}_{\text{ij}}^{o} + {\overrightarrow{v}}_{\text{ij}}^{o}\left( t - t^{o} \right) \right\rbrack \bullet \left\lbrack {\overrightarrow{r}}_{\text{ij}}^{o} + {\overrightarrow{v}}_{\text{ij}}^{o}\left( t - t^{o} \right) \right\rbrack$$                          (Eqn 3-6)
     $$\sigma_{\text{ij}}^{2} = {\overrightarrow{r}}_{\text{ij}}^{o} \bullet {\overrightarrow{r}}_{\text{ij}}^{o} + 2{\overrightarrow{r}}_{\text{ij}}^{o}{\overrightarrow{v}}_{\text{ij}}^{o}\left( t - t^{o} \right) + {\overrightarrow{v}}_{\text{ij}}^{o}{\overrightarrow{v}}_{\text{ij}}^{o}\left( t - t^{o} \right)^{2}$$   (Eqn 3-7)

Here we define a variable $b_{\text{ij}}$ as a dot product,

  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- -----------
     $$b_{\text{ij}} \equiv {\overrightarrow{r}}_{\text{ij}} \bullet {\overrightarrow{v}}_{\text{ij}} = \left| {\overrightarrow{r}}_{\text{ij}} \right|\left| {\overrightarrow{v}}_{\text{ij}} \right|\cos\beta$$   (Eqn 3-8)
  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- -----------

If $b_{\text{ij}} > 0$, the relative position
$\left| {\overrightarrow{r}}_{\text{ij}} \right|$ will be increasing
with time. The two particles will be further away with each other. If
$b_{\text{ij}} = 0$, the relative velocity and relative position vectors
are perpendicular. If $b_{\text{ij}} < 0$, the relative position
$\left| r_{\text{ij}} \right|$ will be decreasing and there may be a
collision. So after removing the superscript "o" on
${\overrightarrow{r}}_{\text{ij}}$ and
${\overrightarrow{v}}_{\text{ij}}$, E can be rewritten as follows,

  -- ----------------------------------------------------------------------------------------------------------------------------------------- -----------
     $$\sigma_{\text{ij}}^{2} = r_{\text{ij}}^{2} + 2b_{\text{ij}}\left( t - t^{o} \right) + v_{\text{ij}}^{2}\left( t - t^{o} \right)^{2}$$   (Eqn 3-9)
  -- ----------------------------------------------------------------------------------------------------------------------------------------- -----------

or

  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------
     $$v_{\text{ij}}^{2}\left( t - t^{o} \right)^{2} + 2b_{\text{ij}}\left( t - t^{o} \right) + \left( r_{\text{ij}}^{2} - \sigma_{\text{ij}}^{2} \right) = 0$$   (Eqn 3-10)
  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------

Solving for time increment by using the quadratic formula, we can have,

  -- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $$t^{c} = t - t^{o} = \frac{- b_{\text{ij}} + s_{\mathrm{\text{time}}}\sqrt{b_{\text{ij}}^{2} - v_{\text{ij}}^{2}\left( r_{\text{ij}}^{2} - \sigma_{\text{ij}}^{2} \right)}}{v_{\text{ij}}^{2}}$$   (Eqn 3-11)
  -- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

where$s_{\mathrm{\text{time}}} = \pm 1$. The discriminant,
$b_{\text{ij}}^{2} - v_{\text{ij}}^{2}\left( r_{\text{ij}}^{2} - \sigma_{\text{ij}}^{2} \right)$,
in E has to be non-negative for a collision to occur. If it is negative,
then the particles will miss each other. The smaller real root
$s_{\mathrm{\text{time}}} = - 1$ is used for the solution to the
quadratic as shown in Figure. The smaller real root occurs when the
particles collide at ${\overrightarrow{r}}_{j}^{s -}$. The location
labeled ${\overrightarrow{r}}_{j}^{s +}$ is the root that will not occur
physically.

*Figure 3-7 Only the smaller real root,*
$s_{\mathrm{\text{time}}} = - 1$ *is used.*

If the particle has a potential shell (See Figure), the equation for the
event time can be modified to,

  -- ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $$t^{c} = t - t^{o} = \frac{- b_{\text{ij}} + s_{\mathrm{\text{time}}}\sqrt{b_{\text{ij}}^{2} - v_{\text{ij}}^{2}\left( r_{\text{ij}}^{2} - d_{\text{ij}}^{2} \right)}}{v_{\text{ij}}^{2}}$$   (Eqn 3-12)
  -- ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

where $d_{\text{ij}}$ is the actual distance of centers for the event.

*Figure 3-8 (a) Approaching particles that will enter the potential
shell; (b) particles in the shell that will have a hard sphere
collision; (c) particles in the shell that will miss a hard sphere
collision, but experience a potential l change event after they pass;
(d) particles moving away that will experience a potential change event
at the shell boundary.*

Calculation of velocity changes upon an event will be a little more
complex. We use prime (') to denote the state after the event. The
change in potential energy for the event is indicated by
$\mathrm{\Delta}\varepsilon = \varepsilon^{'} - \varepsilon$, which may
be zero, positive or negative depending on whether it is a hard sphere
collision, or particles enter or escape from a potential shell.

Total energy is conserved upon an event between particles $i$ and $j$,

  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------
     $${\frac{1}{2}m}_{i}{\overrightarrow{v}}_{i} \bullet {\overrightarrow{v}}_{i} + {\frac{1}{2}m}_{j}{\overrightarrow{v}}_{j} \bullet {\overrightarrow{v}}_{j} + \varepsilon = {\frac{1}{2}m}_{i}{\overrightarrow{v}}_{i}^{'} \bullet {\overrightarrow{v}}_{i}^{'} + {\frac{1}{2}m}_{j}{\overrightarrow{v}}_{j}^{'} \bullet {\overrightarrow{v}}_{j}^{'} + \varepsilon'$$   (Eqn 3-13)
  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------

which can be rearranged,

  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $${\frac{1}{2}m}_{i}\left( {\overrightarrow{v}}_{i} - {\overrightarrow{v}}_{i}^{'} \right) \bullet \left( {\overrightarrow{v}}_{i} + {\overrightarrow{v}}_{i}^{'} \right) = - {\frac{1}{2}m}_{j}\left( {\overrightarrow{v}}_{j} - {\overrightarrow{v}}_{j}^{'} \right) \bullet \left( {\overrightarrow{v}}_{j} + {\overrightarrow{v}}_{j}^{'} \right) + \mathrm{\Delta}\varepsilon$$   (Eqn 3-14)
  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

When the two particles contact, the relative position will be,

  -- -------------------------------------------------------------------------------------------------------- ------------
     $${\overrightarrow{r}}_{\text{ij}}^{c} = {\overrightarrow{r}}_{i}^{c} - {\overrightarrow{r}}_{j}^{c}$$   (Eqn 3-15)
  -- -------------------------------------------------------------------------------------------------------- ------------

The change in momentum for an event is proportional to the collision
vector ${\overrightarrow{r}}_{\text{ij}}^{c}$,

  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $$m_{i}\left( {\overrightarrow{v}}_{i}^{'} - {\overrightarrow{v}}_{i} \right) = - m_{j}\left( {\overrightarrow{v}}_{j}^{'} - {\overrightarrow{v}}_{j} \right) \propto {\overrightarrow{r}}_{\text{ij}}^{c}$$   (Eqn 3-16)
  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

or

  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------
     $$\frac{\left( {\overrightarrow{v}}_{i}^{'} - {\overrightarrow{v}}_{i} \right)}{m_{j}} = - \frac{\left( {\overrightarrow{v}}_{j}^{'} - {\overrightarrow{v}}_{j} \right)}{m_{i}} = \phi{\overrightarrow{r}}_{\text{ij}}^{c}$$   (Eqn 3-17)
  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------

$\phi$ is the parameter we will try to determine. Insert *Eqn 5-14* into
*Eqn 5-11*,

  -- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $$\phi{\overrightarrow{r}}_{\text{ij}}^{c} \bullet \left( {\overrightarrow{v}}_{i} + {\overrightarrow{v}}_{i}^{'} \right) = \phi{\overrightarrow{r}}_{\text{ij}}^{c} \bullet \left( {\overrightarrow{v}}_{j} + {\overrightarrow{v}}_{j}^{'} \right) - \frac{2\mathrm{\Delta}\varepsilon}{\left( m_{i}m_{j} \right)}$$   (Eqn 3-18)
  -- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

Rearrange *Eqn 5-14* also leads to,

  -- ------------------------------------------------------------------------------------------------------------- ------------
     $${\overrightarrow{v}}_{i}^{'} = {\overrightarrow{v}}_{i} + m_{j}\phi{\overrightarrow{r}}_{\text{ij}}^{c}$$   (Eqn 3-19)
  -- ------------------------------------------------------------------------------------------------------------- ------------

and

  -- ------------------------------------------------------------------------------------------------------------- ------------
     $${\overrightarrow{v}}_{j}^{'} = {\overrightarrow{v}}_{j} - m_{i}\phi{\overrightarrow{r}}_{\text{ij}}^{c}$$   (Eqn 3-20)
  -- ------------------------------------------------------------------------------------------------------------- ------------

Plugging *Eqn 5-16* and *Eqn 5-17* into *Eqn 5-15* to eliminate the
primed variables,

  -- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $$\phi{\overrightarrow{r}}_{\text{ij}}^{c} \bullet \left( 2{\overrightarrow{v}}_{i} + m_{j}\phi{\overrightarrow{r}}_{\text{ij}}^{c} \right) = \phi{\overrightarrow{r}}_{\text{ij}}^{c} \bullet \left( 2{\overrightarrow{v}}_{j} - m_{i}\phi{\overrightarrow{r}}_{\text{ij}}^{c} \right) - \frac{2\mathrm{\Delta}\varepsilon}{\left( m_{i}m_{j} \right)}$$   (Eqn 3-21)
  -- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

Define
$d_{\text{ij}}^{2} \equiv {\overrightarrow{r}}_{\text{ij}}^{c} \bullet {\overrightarrow{r}}_{\text{ij}}^{c}$,
${\overrightarrow{v}}_{\text{ij}} \equiv {\overrightarrow{v}}_{i} - {\overrightarrow{v}}_{j}$
and
$b_{\text{ij}}^{c} \equiv {\overrightarrow{r}}_{\text{ij}}^{c} \bullet {\overrightarrow{v}}_{\text{ij}}$,
and rearrange *Eqn 5-18*,

  -- ---------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $$\left( {m_{i} + m}_{j} \right)d_{\text{ij}}^{2}\phi^{2} + 2b_{\text{ij}}^{c}\phi + \frac{2\mathrm{\Delta}\varepsilon}{\left( m_{i}m_{j} \right)} = 0$$   (Eqn 3-22)
  -- ---------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

Applying the quadratic formula,

  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------
     $$\phi = \frac{- 2b_{\text{ij}}^{c} \pm \sqrt{4\left( b_{\text{ij}}^{c} \right)^{2} - 8\left( {m_{i} + m}_{j} \right)d_{\text{ij}}^{2}\mathrm{\Delta}\varepsilon/\left( m_{i}m_{j} \right)}}{2\left( {m_{i} + m}_{j} \right)d_{\text{ij}}^{2}}$$   (Eqn 3-23)
  -- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------

Define the reduced mass, $m_{r} = \frac{m_{i}m_{j}}{m_{i}{+ m}_{j}}$,
then,

  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------
     $$\phi = \frac{- b_{\text{ij}}^{c} + s_{\mathrm{\text{vel}}}\sqrt{\left( b_{\text{ij}}^{c} \right)^{2} - 2d_{\text{ij}}^{2}\mathrm{\Delta}\varepsilon/m_{r}}}{\left( {m_{i} + m}_{j} \right)d_{\text{ij}}^{2}}$$   (Eqn 3-24)
  -- ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ------------

where the sign of $s_{\mathrm{\text{vel}}}$ depends on the context of
the event (See Figure). For an event where
$\mathrm{\Delta}\varepsilon = 0$, $s_{\mathrm{\text{vel}}} = - 1$ is the
only reasonable solution. For a hard sphere collision event, the
distance is the inter-particle radius,
$d_{\text{ij}} = \sigma_{\text{ij}}$ and the relation becomes,

  -- ---------------------------------------------------------------------------------------------- ------------
     $$\phi = \frac{- 2b_{\text{ij}}^{c}}{\left( {m_{i} + m}_{j} \right)\sigma_{\text{ij}}^{2}}$$   (Eqn 3-25)
  -- ---------------------------------------------------------------------------------------------- ------------

Based on the context of the event, the DMD calculation will follow the
flowsheet below,

![C:\\Users\\Size\\Desktop\\22.png](media/image8.png){width="6.0in"
height="4.116029090113736in"}

*Figure 3-9 Flowsheet for DMD with wells.^9^*

### **Energy Minimization**

Theoretically, the sDMD can import any PDB file of proteins to create
the simulation system. In the PDF file, the original coordinates of
atoms may violate the ranges of bond lengths or the minimum distances
between pairs of atoms in the sDMD\'s library. To eliminate these
conflicts, the sDMD uses steepest descent algorithm.

In DMD, the potential function is a step function, then the first-order
derivative would not be possible. However, by computing the relative
positions of a pair of atoms, we could obtain the vectors that can guide
the atoms to move closer or further; based on their distance, then, we
can find at which step the potential currently locates, as well as to
which direction the potential will decrease. From these information, we
can decide the next move of an atom to approach the local minimum of
potential energy.

In practice, to achieve a fast energy minimization, the sDMD will only
eliminate the bond length violations and atom diameter overlapping. In
other words, the program will only vibrate an atom if it has bond that
was too short or too long and if it overlaps any surrounding atoms. User
can follow a short simulation under a low temperature to further relax
the system.

### **Confined Walls**

The sDMD is able to simulate a system in various types of confined
environments, such as parallel pores, cylindrical tubes, and spherical
cavities. Currently the confined walls can be only smooth surfaces,
which means there are no explicit atoms on the walls. The wall surface
is like a multi-layer shell, and each layer represents a potential step.
Or we can image the wall is made of extremely high density of atoms, so
high that the other atoms interact with it will be like interacting with
a smooth surface.

To simulate the parallel pores is very straightforward: once the
positions of the walls are set, the event time in the nearest future
when an atom interact with the wall can be easily determined, and the
after-interaction velocity can be computed directly, since it is just a
one-dimension problem. However, for a cylindrical tube or a spherical
cavity, the situation will become a little bit more complex. Due to the
curvature of surface and the implicit atoms, it would be difficult to
determine the exact direction and position on the wall at which the atom
will interact at the nearest event, which both are required to compute
the interaction potential and after-interaction velocity change. To deal
with this issue, the sDMD uses an indirect way: it is easy, instead, to
compute the instant distance of an atom from the center of either a
cylinder tube or a spherical cavity, L. Then the distance between the
atom and the wall will be D = R -- L, where R is the diameter of the
tube or cavity. Thus, the potential range in the potential step can be
determined by D. Next, instead of doing a lot of trigonometric
calculations or trying to find a pseudo-atom on the wall to model the
interaction with the target atom, the sDMD consider this as an event of
interaction between the target atom and a huge pseudo-atom at the center
of the tube or cavity. In other words, any interactions between an atom
and the wall will be like an atom inside another huge atom interacts
internally with the shell of the huge atom. The interaction types will
be switched. For example, the well-capture event corresponding to an
atom and the center pseudo-atom will become a well-escape event
corresponding the atom and the wall, and similarly, the hard-sphere
collision event will become a well-bounce, etc. Thus, the sDMD could use
the same algorithm to compute the interaction between a pair of atoms
and between atom and wall.

![](media/image9.jpg){width="2.5in" height="2.5in"}

*Figure 3-10 Schematic diagram of wall interaction in a cylindrical tube
or a spherical cavity. The big solid circle is the wall boundary, which
is infinitely repulsive; the big dashed circle is the potential shell of
wall. The small black dashed circle is the original position of the
target atom; the small black solid circle is the future position of the
target atom after time t. The small read circle is the pseudo-atom at
the center. The black line is the diameter of the wall, R; the red line
is the distance between the target atom and the pseudo-atom, L; the blue
line is the distance between the target atom and the wall, D = R - L.*

The sDMD also allows users to use the obstruction wall. The obstruction
wall is an infinite wall parallel to the xy-, xz- or yz-planes. Multiple
obstruction wall can be inserted into the system box at any positions.
This can be used to construct close-end cylindrical tubes, cuboid
confined environments, and tunnel and reservoirs (combined with tunnel
function), see Figure 3-11, or divide the system box into different
regions.

![C:\\Users\\Size\\AppData\\Local\\Microsoft\\Windows\\INetCache\\Content.Word\\Flow.png](media/image10.png){width="3.5208333333333335in"
height="2.0in"}

*Figure 3-11 Combination of obstruction walls and cylindrical tube
confinement can create a tube between two reservoirs, which can be used
to model protein translocations. The solid rectangle represents the
boundary of the system box; the black dashed rectangles represent the
obstruction walls; the red dashed rectangle represents the cylindrical
tube.*

### **Replica Exchange Molecular Dynamics Simulation (REMD)**

The sDMD allows users to use REMD. The results can be analyzed by
T-WHAM^10-11^ in the analysis executable.

REMD can be used to efficiently explore the potential energy landscape
of molecular systems. The ruggedness and the slope toward the energy
minimum in the landscape govern sampling efficiency at a given
temperature. Although escape out of local minima is accelerated at
higher temperatures, the free energy landscape is altered due to large
entropic contributions. To efficiently overcome energy barriers while
maintaining conformational sampling corresponding to a relevant free
energy surface, the sDMD employs the replica exchange sampling
scheme^8,\ 12-13^.

In REMD, multiple simulations/replicas of the same system are performed
in parallel at different temperatures. Individual simulations are
coupled through Monte Carlo-based exchanges of simulation temperatures
between replicas at periodic time intervals. Temperatures are exchanged
between two replicas, i and j, maintained at temperatures Ti and Tj and
with energies Ei and Ej, according to the canonical Metropolis criterion
with the exchange probability $P$, where $P = 1$ if
$\mathrm{\Delta} = \left( \frac{1}{k_{B}T_{i}} - \frac{1}{k_{B}T_{j}} \right)\left( E_{j} - E_{i} \right) \leq 0$,
and $P = e^{- \mathrm{\Delta}}$, if $\mathrm{\Delta} > 0$. To achieve
this, the sDMD uses a Socket server to send and receive the temperature
and energy data between pairs of replicas at the preset time intervals.
The sDMD would suggest the exchange rate to be every 1 x 10^3^ time
units.

1\. Voegler Smith, A.; Hall, C. K., alpha-helix formation: discontinuous
molecular dynamics on an intermediate-resolution protein model.
*Proteins* **2001,** *44* (3), 344-60.2. Voegler Smith, A.; Hall, C. K.,
Bridging the gap between homopolymer and protein models: A discontinuous
molecular dynamics study. *The Journal of Chemical Physics* **2000,**
*113* (20), 9331-9342.3. Smith, A. V.; Hall, C. K., Assembly of a
tetrameric alpha-helical bundle: computer simulations on an
intermediate-resolution protein model. *Proteins* **2001,** *44* (3),
376-91.4. Smith, A. V.; Hall, C. K., Protein refolding versus
aggregation: computer simulations on an intermediate-resolution protein
model. *J Mol Biol* **2001,** *312* (1), 187-202.5. Buldyrev, S. V.,
Application of Discrete Molecular Dynamics to Protein Folding and
Aggregation. In *Aspects of Physical Biology: Biological Water, Protein
Solutions, Transport and Replication*, Franzese, G.; Rubi, M., Eds.
Springer Berlin Heidelberg: Berlin, Heidelberg, 2008; pp 97-131.6. S.
Takada, Z. L.-S., P. G. Wolynes, Folding dynamics with nonadditive
forces: A simulation study of a designed helical protein and a random
heteropolymer. *J. Chem. Phys.* **1999,** *110* (23), 11616-11629.7.
Andersen, H. C., Molecular dynamics simulations at constant pressure
and/or temperature. *The Journal of Chemical Physics* **1980,** *72*
(4), 2384-2393.8. Ding, F.; Tsao, D.; Nie, H.; Dokholyan, N. V., Ab
initio folding of proteins with all-atom discrete molecular dynamics.
*Structure* **2008,** *16* (7), 1010-8.9. Elliott, J. R.; Lira, C. T.,
Supplement. In *Introductory chemical engineering thermodynamics*, 2nd
ed.; Prentice Hall: Upper Saddle River, NJ, 2012.10. Kumar, S.;
Rosenberg, J. M.; Bouzida, D.; Swendsen, R. H.; Kollman, P. A., THE
weighted histogram analysis method for free-energy calculations on
biomolecules. I. The method. *Journal of Computational Chemistry*
**1992,** *13* (8), 1011-1021.11. Gallicchio, E.; Andrec, M.; Felts, A.
K.; Levy, R. M., Temperature weighted histogram analysis method, replica
exchange, and transition paths. *J Phys Chem B* **2005,** *109* (14),
6722-31.12. Zhou, R.; Berne, B. J.; Germain, R., The free energy
landscape for beta hairpin folding in explicit water. *Proc Natl Acad
Sci U S A* **2001,** *98* (26), 14931-6.13. Okamoto, Y.,
Generalized-ensemble algorithms: enhanced sampling techniques for Monte
Carlo and molecular dynamics simulations. *J Mol Graph Model* **2004,**
*22* (5), 425-39.

**Program Structure**
=====================

This section will briefly describe the structure of the program.

struct AtomStr{} is the most important data structure in the sDMD. It
stores almost all the information of an individual atom.

"property" stores all the static information, like mass and atom number;
"dynamic" stores all the dynamic information, like coordinate and
velocity; "eventList" stores the position of this atom in the search
binary tree, as well as its parent and children.

The sDMD has been optimized by the Pthreads Library for the future
development of parallel computing. The program will establish
single/multiple threads and each of the thread can perform the
calculation by its own. Thus, the data in the sDMD have two layers, one
is the raw data and the other is the data copy in each thread. Figure
4-1 shows a diagram of a single-core simulation.

![](media/image11.jpg){width="3.996583552055993in" height="2.5in"}

*Figure 4-1 Single-core simulation process*

At the beginning of a simulation, the thread will be assigned a target
atom whose event is in the nearest future. If the target atom's event is
an interaction event, the thread will be also assigned a partner atom
with which the target atom will interact. The thread will make a copy of
the data of these atoms and do the following calculations by only using
this copy. After the event being processed, the thread will further
predict the next possible event for this target atom. During this whole
process, the raw data will not be affected. Before committing the
processed data and the predicted information to the raw data, the
program will check the hazards. The hazard here means any situations
that would make the previously predicted event outdated. In other words,
if there exists any hazard, the processed event by the thread will be
invalid. If this is the case, the thread will re-calculate the next
event for the target atom or, in some cases, would be assigned by
another target atom and do the same calculations. Only if there exists
no hazard, the calculated results will be committed to the raw data.

The sDMD will keep looping the above processes until the assigned
simulation time is reached.

5.  **Program Usage**
    =================

    3.  **Installation**
        ----------------

The package includes three folders: simulation source code folder,
"src", analysis source code folder, "analysis", and an "input" folder.
Each source code folder has its dedicated Make file, "Makefile". The
"input" folder contains the configuration file, "parameter.txt", and the
force filed library folder, "Library".

To compile the simulation source code, use the following command,

*\$make all*

It will generate three executables, sDMD, sServer, and sREMD. sDMD is
the core program to run the simulation; sREMD is the program to launch
REMD simulation, only; sServer is an assistant program to create a
socket server for REMD simulation, and will be executed by sREMD
automatically.

If users only need the core program to run regular DMD simulations,
which only the sDMD executable will be required, simply use the command,

*\$make sDMD*

To compile the analysis source code, use the following command,

*\$make*

It will generate only one executable, analysis.

4.  **Run Parameters**
    ------------------

    7.  ### **Non-REMD simulation**

sDMD will need a coordinate file, \*.gro or \*.pdb, the configuration
file, parameter.txt, and the force field library folder, "Library", to
start a simulation.

Below is an example of the configuration file,

  continue?      If users would like to perform a new simulation, use "new"; if users would like to restart a previous run, which may have not finished yet, or would like to extend the current simulation, use "continue". Both process will back up the old files. The backup files will have a suffix representing the backing up time and date.
  -------------- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Time           Simulation time.
  T(reduced)     Reduced system temperature, T\* = T / 502 K. If set to "no", the program will keep using the temperature set at the last time.
  Dump\_Rate     How often the data will dump. Here is each 10 time units.
  CutOff\_R(A)   The cutoff radius of VDW interaction. Unit is A.
  Force\_Field   The force field employed in the simulation. Currently only support "Ding"^1^.
  Thermostat     The thermostat algorithm. Currently only support "Andersen"^2^.
  PDMD\_Method   Parallel algorithm. Currently only support "1", which means only one CPU core will be used.
  ThreadNo       Number of cores employed to perform the parallel computing. Currently only support "1".
  WallExist?     If users would like to perform a simulation under bulk conditions, use "no". If users would like to perform a simulation in confined environments, use "smooth".
  WallType       Shape of the wall. Support "parallel", "cylinder", and "sphere".
  WallDyn        If the confined walls are fixed, use "no". If the confined walls are dynamic, use "compress", "extend", or "repulse".
  Flow(x,y,z)    If there will be constant velocities applied on the atoms. Each on x-, y-, and z-axis.
  Obstruct?      If there exists any obstruction walls in the system. If no, use "no". Otherwise, the first number represents how many of the obstruction walls will be inserted into the system, following the position set of each obstruction wall on x-, y-, and z-axis, separating by commas.
  Tunnel?        If there will be a tunnel in the system. This is only used for the translocation simulation. If no, use "no". Otherwise, the three numbers represent, respectively, the start coordinate, the end coordinate and the diameter of the tunnel. These coordinates are all on x-axis and the unit is A.

To run the executable sDMD, use the following command,

\$./sDMD \[FLAGS\]

Use the flag -h or -help to show the supported FLAGs,

+-----------------------------------+-----------------------------------+
| -i                                | directory of the parameter folder |
+===================================+===================================+
| -o                                | directory of the output folder    |
+-----------------------------------+-----------------------------------+
| -f                                | (optional) follow the file name   |
|                                   | of the input coordinate file      |
+-----------------------------------+-----------------------------------+
| -si                               | (optional) follow the file name   |
|                                   | of the input saved data for       |
|                                   | continuity                        |
+-----------------------------------+-----------------------------------+
| -so                               | (optional) follow the file name   |
|                                   | of the output saved data for      |
|                                   | continuity                        |
+-----------------------------------+-----------------------------------+
| -trj                              | (optional) follow the file name   |
|                                   | of the output trajectory file     |
+-----------------------------------+-----------------------------------+
| -cnt                              | (optional) follow the file name   |
|                                   | of the output connection map      |
+-----------------------------------+-----------------------------------+
| -log                              | (optional) follow the file name   |
|                                   | of the output Log file            |
+-----------------------------------+-----------------------------------+
| -sys                              | (optional) follow the file name   |
|                                   | of the output SysInfo file        |
+-----------------------------------+-----------------------------------+
| -pot                              | (optional) output potential       |
|                                   | energy, follow the file name,     |
|                                   | otherwise will use the default    |
|                                   | name                              |
+-----------------------------------+-----------------------------------+
| -kin                              | (optional) output kinetic energy, |
|                                   | follow the file name, otherwise   |
|                                   | will use the default name         |
+-----------------------------------+-----------------------------------+
| -tem                              | (optional) output temperature,    |
|                                   | follow the file name, otherwise   |
|                                   | will use the default name         |
+-----------------------------------+-----------------------------------+
| -HBn                              | (optional) output hydrogen bond   |
|                                   | number, follow the file name,     |
|                                   | otherwise will use the default    |
|                                   | name                              |
+-----------------------------------+-----------------------------------+
| -xyz                              | (optional) output xyz trajectory, |
|                                   | follow the file name, otherwise   |
|                                   | will use the default name         |
+-----------------------------------+-----------------------------------+
| -pdb                              | (optional) output pdb trajectory, |
|                                   | follow the file name, otherwise   |
|                                   | will use the default name         |
+-----------------------------------+-----------------------------------+
| -box                              | (required only if the coordinate  |
|                                   | file is a .PDB file) the          |
|                                   | simulation box dimensions, x, y,  |
|                                   | z                                 |
+-----------------------------------+-----------------------------------+
| -Wsz                              | (required only if WallDyn is      |
|                                   | assigned) follow the maximum      |
|                                   | changing size of the wall,        |
|                                   | default 5 A                       |
+-----------------------------------+-----------------------------------+
| -Wrt                              | (required only if WallDyn is      |
|                                   | assigned) follow the total time   |
|                                   | needed for the size change,       |
|                                   | default 200                       |
+-----------------------------------+-----------------------------------+
| -REMD                             | only use during REMD, follow the  |
|                                   | server name,                      |
|                                   |                                   |
|                                   | the port number,                  |
|                                   |                                   |
|                                   | the temperature,                  |
|                                   |                                   |
|                                   | the exchange rate, and            |
|                                   |                                   |
|                                   | the suffix name of saving files   |
|                                   |                                   |
|                                   | (this flag will be set            |
|                                   | automatically by REMD executable) |
+-----------------------------------+-----------------------------------+

### **REMD Simulation**

Besides the required files and commands mentioned above, to run a REMD
simulation, the sDMD needs another configuration file, "REMDConf.txt",
and the other two executables, server and REMD. The configuration file
must be in the "input" folder and the two executables must be in the
same folder of executable sDMD. The "T(reduced)" item should be set to
"no".

Below is an example of the configuration file (the first line is just a
comment),

  ReplicaNum     The number of replicas.
  -------------- ---------------------------------------------------------------------
  Temperatures   The temperature of each replica.
  SocketPort     The port number of the socket server.
  ExchangeRate   How often the temperature will be exchanged. Unit is the time unit.

To run the executable REMD, use the following command,

\$./REMD \[distribute T or not\] -f \[configuration file\] -args \[args
of executable without flag -REMD\]

  -nodist        \[yes\] Default is yes, if absent. Otherwise, the program will NOT distribute the preset temperature to each replica. Use this flag if users would like to extend/restart the simulation.
  -------------- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Temperatures   The temperature of each replica.
  SocketPort     The port number of the socket server.
  ExchangeRate   How often the temperature will be exchanged. Unit is the time unit. Suggest to be 1000.

The flags of executable server and the -REMD flag of executable sDMD
will be set by REMD automatically.

### **Analysis**

To perform an analysis, use the following command,

\$./analysis \[FLAGS\]

Use the flag -h or -help to show the supported FLAGs,

  -path     exact path of data folder
  --------- ------------------------------------------------------------------------
  -trj      trajectory input file
  -cnt      connect map input file
  -sys      (optional) system info input file
  -log      (optional) log input file
  -sum      (optional) ignore the old time record in .log file
  -rPBC     remove PBC, call at the beginning. require -trj input file
  -HB       analyze HB info. require -cnt input file
  -En       analyze energy info
  -Ag       analyze aggregation info
  -RG       calculate RG (not support yet)
  -MSD      calculate MSD (not support yet)
  -Ramach   (require -HB) plot Ramachandran plots for each amino acid
  -ConMap   (require -HB) plot contact maps for atoms and amino acids
  -REMD     analyze REMD data. follow the flags below:
  -reNo     number of replicas
  -bin      (optional) number of bins in histogram. default 10
  -aHB      (optional) use alpha HB number as the reaction coordinate
  -bHB      (optional) use beta HB number as the reaction coordinate
  -tHB      (optional) use total HB number as the reaction coordinate
  -maxE     (optional) max value of potential energy (most negative). default -150
  -minE     (optional) min value of potential energy (least negative). default 0
  -temp     (optional) target temperature. default 300
  -st       (optional) start frame of analyzing. default 0
  -et       (optional) end frame of analyzing. default -1 (the final frame)

The analysis program is also highly modularized and optimized. It will
be easy to add more analysis functions onto it.

6.  **Format**
    ==========

    5.  **Coordinate File**
        -------------------

The sDMD supports .pdb and .gro files as the coordinate input file. But
it recommends to use a .gro file as the coordinate file, since the .pdb
file has too many alternative formats.

When the sDMD imports the coordinate file, it will first read the
original .pdb or .gro file. During this process, the program will remove
the atoms those are not in the model library, such as the light
hydrogens and one of the carboxylic oxygen on the tail amino acid. At
the same time, the program will extract the atom information and assign
them onto the variables about protein/peptide, amino acid and atom.
After read through the whole coordinate file, some data will be dumped
into a .xyz file. The program will then read the .xyz file in order to
finish the final assignments for some other variables. The reason of
such a process is that some of the variables can be correctly assigned
only after reading through and modifying the original coordinate file.

The output trajectory file will be of .gro format. It consists of
multiple continuous snapshots of the atoms' coordinates.

### **GRO File**

Pay attention that the coordinate unit of .gro file is nm. The length
unit in the sDMD, however, is A. The program will convert the unit
automatically.

The following is an example of the .gro file,

The first line is a comment, can be a brief description of the file.
Only one line of comment will be allowed.

The second line has only one number represents the total atom number in
the system. In this example it has 375 atoms. To save pages it does not
show them all.

The line after the second contains the information of each atom. Use the
line below as an example,

2ALA N 9 1.256 2.353 3.680

The first term describes at which amino acid in a protein/peptide the
atom locates. Here is the 2^nd^ amino acid whose name is "ALA". Pay
attention that the amino acid number here is not continuous between
molecules, which means that if there are multiple molecules in the
system, this number will start from 1 for each molecule. This will help
the program to judge how many molecules there are in the system. The
second term is the atom's name. The third term is the number of the atom
in the whole system. This number is continuous in the whole system. The
following three terms are the atom's coordinates on x-, y- and z-axis,
respectively. If this is a output .gro file from the simulation, there
will be three more terms following the coordinates. They are the atom's
velocities on x-, y- and z-axis, respectively.

The three numbers on the last line of the .gro file are the sizes of the
system box, on x-, y- and z-axis, respectively.

### **PDB File**

The following is an example of the supported .pdb file,

The .pdb file can have several lines of comments/descriptions at the
head. The useful information is on the line with the keyword "ATOM" in
front. Use the following line as an example,

ATOM 9 N ALA A 2 15.270 25.150 38.090 1.00 0.00 N

The first term is the keyword. The second term is the number of the atom
in the whole system, which is continuous. The third term is the atom's
name. The fourth term is the amino acid name at which the atom locates.
The fifth term is just an arbitrary symbol for the protein molecule. The
atoms in the same molecule will have the same symbol. Different
molecules can share the same symbol, but not the continuous two. The
sixth term is the amino acid's number in the protein. This number is not
continuous between different molecules. The following three terms are
the atom's coordinates on x-, y- and z-axis, respectively. The sDMD will
ignore all the terms behind.

To import from the coordinate file of a .pdb format, users need to
specify the box size by using -box flag while executing sDMD.

**Force Field File**
--------------------

The force field consists of two types of files. One contains the
parameters of amino acid and the other contains the potential table.

### **Parameters of amino acid**

There are 21 parameter files for each of the amino acid (PRO is not
supported yet). The file's name indicates which amino acid it
represents. Here use ALA.txt as an example,

The first line is the name of the amino acid. The lines starting with
the keyword "ATOM" list all the atoms in this amino acid. Pay attention
that all the amino acids, except for "ALA" and "GLY", share the same
backbone characteristics, so for those amino acids the "ATOM" lines only
list their own dedicated atoms and the backbone data will be read from
"ALA.txt" directly. Use the following line as an example,

ATOM CA CA1 12.0 N N

The second term is the atom's name. The third term is the atom's type.
(All types of atoms can be find in "InteractionPotentialTable.txt" file
in the library.) The fourth term is the atom's mass. Pay attention that
the mass of hydrogen is set as 12.0, as the author suggests.^1^ In some
cases this will have the hydrogen bond harder to form. In the sDMD, it
changes the mass back to 1.0. The user can disable it by changing the
source code, function "ScanAA" in "initialization.c". The last two terms
are left for the further development. They can be assigned by any values
which can give the atom extra properties.

The lines starting with the keyword "BOND" list all the bonds in this
amino acid. Again, all the amino acids except "ALA" and "GLY" only list
their own dedicated bonds. Use the following line as an example,

BOND N H 1.000 0.020

This line means there is a covalent bond between atom N and atom H. The
bond length is 1.000 A and can vibrate between
$1.0 \times \left( 1 - 0.02 \right)$ and
$1.0 \times \left( 1 + 0.02 \right)$.

The format of the lines starting with the keywords "ANGLE" and "INTER"
are similar with "BOND" lines. They are however auxiliary bonds. The "+"
symbol in front of the atom name means this atom locates in the next
amino acid.

The line starting with the keyword "CONSTR" list all the dihedral
constraints in the amino acid. Use the following line as an example.

CONSTR2 C +C 2.700 2.846 1.500 2.910 0.500 4.000

The number after the keyword "CONSTR" represents the number of steps
there will be in the potential function. The second and the third terms
are the head and tail atoms in the dihedral angle. Starting from the
fourth term are the values construct the potential function. The first
and the last values represent the minimum and the maximum distance
between the two atoms; the values in the middle are presented as pairs:
the first values are the discontinuity positions in the potential
function and the second values are the energy change crossing the steps.
In this example, the maximum distance is 4.000 A; when the two atoms
move closer than 2.910 A, the potential energy will decrease 0.500 kcal
/ mol, or become less negative by 0.500 kcal / mol; when they move even
closer than 2.846 kcal / mol, the potential energy will further decrease
1.500 kcal / mol. The distance between the two atoms cannot be smaller
than 2.700 A.

Currently the sDMD only allows to use one type of atom to construct the
wall. The atom used is listed in the file "Wall.txt". The default is
"CA". Users can change it at will.

### **Potential Table**

There are two types of potential table: one is for VDW interaction,
called "InteractionPotentialTable.txt", and the other is for hydrogen
bond interaction, called "HBPotentialTable.txt". They share the same
format and are similar to that used in "CONSTR" line. Any lines starting
with "\#" symbol are the comments. Use the following line from
"InteractionPotentialTable.txt" as an example,

HB CA 2.000000 2.50000 0.20000 4.50000 0.02000 6.50000 0.01000

The first two terms represent the pair of atoms that holds the following
potential function. The only difference from the "CONSTR" line is that
there is no maximum distance, which means the number of values starting
from the third term is always odd, while in the "CONSTR" line is even.

**Functions**
=============

This section will briefly describe the core functions in the main source
code files. Most functions can be identified by their names.

**Simulation Source Code**
--------------------------

#### DMD.c

This file contains the main function.

#### DataSave.c

The functions in this file manage the naming, backing up, creating,
initializing, saving and closing of simulation files. All the saving
processes will call the function SaveData().

#### Event.c

This file contains all the functions used to calculate the processes of
different events. The DoEvent() manages which event function will be
called. The events like collision, bonding, formation and breaking of
hydrogen bonds all will call InteractionEvent(). Hydrogen bonding also
has its own event function, HBEvent(). Other events all have their own
event functions: ThermostatEvent() for thermostat event,
PBCandCrossCellEvent() for period boundary crossing and cell crossing
event, WallEvent() for regular wall interaction event, ObstEvent() for
obstruction wall interaction event and TunnelEvent() for wall
interaction event in tunnel.

The function LinkList() manages the linked list between cells. To
increase the calculation speed, the sDMD divides a system box into many
cells and keeps tracking the cell number at which an atom locates. Thus,
by using the linked list between the atom number and cell number, it
will be easy to find which atoms are in the target cell. Combining the
linked list and the cut off radius, it will be easy and fast to scan the
27 cells surrounding the target atom by which its event can be
scheduled.

#### Initialization.c

The functions in this file are charged to read and import all the
configurations, parameters and coordinates. They will also allocate the
memory for most of the variables used during the simulation.

The functions are designed based on the formats of the force field.

#### List.c

This file contains the functions to manage the dynamic linked list.

#### Models.c

The functions in this file are used to convert the models' names (of
atoms and amino acids) into numbers.

#### REMD.c

The functions in this file are only used to perform the calculations for
just one replica in a REMD simulation. REMD() initializes the variables.
REMDRun() does the main jobs. ReplicaExchange() manages the data change
and client\_open() initialize a socket client for this replica.

There are two other files to create replicas and server, respectively.

#### SGThread.c

This is the core file to perform a single-core non-REMD simulation.
SingleThread() initializes the variables and SingleThreadRun() does the
main jobs.

#### ThreadProcess.c

The functions of the main procedures to perform a DMD simulation are in
this file. FirstRun() is only used when starting a fresh new simulation.
It will schedule the events for all the atoms and create a search binary
tree to manage the events. ProcessEvent() is used by the thread to
pre-calculate the event, by the function DoEvent(), and pre-predict the
future event, by the function Predict(). These two procedures will be
done only if there exists no hazard for the current event, checked by
the function HazardCheck(). Otherwise, the current event is invalid and
re-prediction will be performed. CommitEvent() is used to commit the
updated atom data to the raw data. AssignJob() is used to assign a
target atom to the thread and AssignThread() is used to update the data
of the thread for the assigned target atom.

#### TimePrediction.c

This file contains all the functions used to predict the event time.
Different events have their own time prediction functions. The function
name indicates the event type.

#### ToolFunctions.c

This file contains all the helper functions. The function FindPair() is
used to find the right potential step at which the target pair of atoms
locates. The Scheduling\*() functions help to manage the search binary
tree of events.

#### sREMD.c

This file is only used to perform the REMD simulation. It helps to
establish multiple replicas and import arguments and flags into the
executables sDMD and sServer.

#### sServer.c

This file is only used in the REMD simulation. It helps to create and
maintain a socket server for the data exchanging between replicas.

**Analysis Source Code**
------------------------

#### Analysis.c

This file contains the main function.

#### FileManage.c

The functions in this file help to assign names to files.

#### SystemInformation.c

The functions in this file are used to read the system information data
from the .dat file and .log file output from the simulation.

#### Cluster.c

The functions in this file are used to calculate the cluster evolution
during aggregation.

#### Energy.c

Th functions in this file are used to calculate the total kinetic
energy, CalKeEnergy(), and potential energy, CalPoEnergy(), of the
system, as well as the potential energy between the molecules and
confinement walls, CalWlEnergy(). CalTemp() is used to calculate the
system temperature, and FindPair() and RightPair() are used to find the
right potential step at which a pair of atoms locates.

#### HBRamach.c

The functions in this file are used to calculate the numbers of
different types of hydrogen bonds. RamachandranPlot() is used to compute
the Ramachandran plot. If the flag "-ConMap" is set, HBInfo() will also
produce the contact maps for both pairs of amino acids and pairs of
atoms.

#### PBCAdjust.c

The functions in this file are used to remove the period boundary
condition in the trajectory file.

#### REMD.c

The functions in this file are used to perform WHAM^3-4^ analysis for
the results from REMD simulations.

#### Tools.c

This file contains all the helper functions.

1\. Ding, F.; Tsao, D.; Nie, H.; Dokholyan, N. V., Ab initio folding of
proteins with all-atom discrete molecular dynamics. *Structure*
**2008,** *16* (7), 1010-8.2. Andersen, H. C., Molecular dynamics
simulations at constant pressure and/or temperature. *The Journal of
Chemical Physics* **1980,** *72* (4), 2384-2393.3. Kumar, S.; Rosenberg,
J. M.; Bouzida, D.; Swendsen, R. H.; Kollman, P. A., THE weighted
histogram analysis method for free-energy calculations on biomolecules.
I. The method. *Journal of Computational Chemistry* **1992,** *13* (8),
1011-1021.4. Gallicchio, E.; Andrec, M.; Felts, A. K.; Levy, R. M.,
Temperature weighted histogram analysis method, replica exchange, and
transition paths. *J Phys Chem B* **2005,** *109* (14), 6722-31.
