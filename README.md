# UMAT for soft tissues with permanent deformation

This constitutive model was developed to simulate the behavior of **soft tissues under cyclic loading**. 

It was seen experimentally that soft tissues can suffer damage under very few loading cycles, a damage mechanism called **ultra-low cycle fatigue** [1][2]

Experiments in sheep pelvic floor muscle shown microdamage in the muscle fibers and visible **permanent deformation (~0.5 permanent strain)** after just 60 cycles at 60% ultimate tensile displacement.

The code is written in FORTRAN and the subroutine can be used with the finite element software ABAQUS. More details regarding the implementation and specific equations can be found in my paper [3].

[1]
[2]
[3]

#Example 2 - Rectangular Specimen under biaxial cyclic loading

![biaxial](https://user-images.githubusercontent.com/95075305/170662306-0929fe3b-2ff6-4e65-baf5-0bd01e6e066b.png)


- Only 1/8 of the specimen was modeled and simetry boundary conditions were applied
- The specimen was subjected to 20 loading cycles where each cycle consists of:
    - stretch to 1.2
    - unstretch to 1.1
- In the last cycle, the specimen is allowed to unload without an imposed displacement so that we can see the permanent deformation.


![biaxial_compressed](https://user-images.githubusercontent.com/95075305/170662729-2d1d3df6-45d1-4efc-9622-0accaad083c9.gif)
