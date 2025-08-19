Scripts used for both generating MD simulations and also to submit jobs in ACENET. 

## Making a simulation using GROMACS
For the sake of not forgetting in the future the general (and also some specifics) steps of using GROMACS for generating a Molecular Dynamics simulation, here I label the steps to follow. My knowledge in this is a combination of making the [Lyzosyme Tutorial](http://www.mdtutorials.com/gmx/lysozyme/index.html), and getting my hands dirty trying to replicate the system of Liu et al.

Additionally to this tutorial (pretty much the basics of GROMACS are laid there), here are some resources that helped me understand better the functioning of GROMACS (and also some of the general concepts in Molecular Dynamics):
 - [Walkthrough of the Lyzosyme Tutorial](https://www.youtube.com/watch?v=rYZ1p5lXNyc&pp=ygUHZ3JvbWFjcw%3D%3D): In case you struggle with the online one (is pretty simple to be honest).
 - [Webinar of GROMACS](https://www.youtube.com/watch?v=MWafKFVgFTU&t=2073s&pp=ygUHZ3JvbWFjcw%3D%3D): Understand the basics of GROMACS.
 - [Presentation on GROMACS](https://www.youtube.com/watch?v=KEfMuHMTBQU&t=1617s&pp=ygUHZ3JvbWFjcw%3D%3D).
 - [GROMACS DOC](https://manual.gromacs.org/current/user-guide/getting-started.html).
 - If you search `gmx` and then the command youre using (for example `gmx mdrun`), you find the documentation of that command, which is pretty useful for understanding the workflow also.
 - [.mdp Parameters DOC](https://manual.gromacs.org/current/user-guide/mdp-options.html): Essential and to understand the inputs of your MD.

All credits to them, this is mostly for myself (so that I can also come back to my notes) but I hope that this can also be useful for someone else starting from zero as I did.

### Setting Up a System

Before doing any simulation, you first need to prepare your system. That is, what you want to simulate, and how the system is defined. Im gonna do a walkthorugh of the collagen fibril I am simulating, as an example, but obviously most of the specific parameters/tweaks need to be made accordingly to your system needs/desires.

#### 1. Making the accordingly GROMACS files

You will start with a PDB file, which has the information of your protein (for instance, for the collagen fibril I am using, I get the file from [ColBuilder](https://colbuilder.h-its.org/), which gives me a PDB file). This PDB file by itself cannot be processed through GROMACS, so the first step is to get the accordingly GROMACS files for your simulation. For this, you use the `gmx pdb2gmx` command (its sort of self explanatory):

```
gmx pdb2gmx -f ticf.pdb -o ticf_p.gro -water tip3p -ff amber14sb
```

Here, you have your input file, output file, the selection of the water model you want to use, and the force field you want to use. For understanding better each an every command step, highly recommend looking at the `gmx ...` documentation online. In this case, we get a geometry file `.gro`, and additionally to the PDB file, we input 2 important things:
 - **The water model**: This is how the water molecules in your system are gonna be "managed"/simulated, for instance that the bonds are flexible, the bonds parameters, as 3 particles or less/more, etc. This parameter is choosen accordingly to experimental data, that is, how well certain water model goes with your protein, from comparison with experimental data. You should look this up, for example for collagen water TIP3P is recommended.
 - **The force field**: Another important parameters that *defines* your simulation. The force field file (`.ff`, is actually a folder) contains all the information regarding the parameters of interaction that define the bond interaction, and the intermolecular interactions, of each particle. 

 

