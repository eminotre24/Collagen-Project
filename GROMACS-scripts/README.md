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
