# UMAT for soft tissues with permanent deformation

This constitutive model was developed to simulate the behavior of **soft tissues under cyclic loading**. 

It was seen experimentally that soft tissues can suffer damage under very few loading cycles, a damage mechanism called **ultra-low cycle fatigue** [1][2]

Experiments in sheep pelvic floor muscle shown microdamage in the muscle fibers and visible **permanent deformation (~0.5 permanent strain)** after just 60 cycles at 60% ultimate tensile displacement.

The code is written in FORTRAN and the subroutine can be used with the finite element software ABAQUS. More details regarding the implementation and specific equations can be found in my paper [3].

###### [1] Chen, J., Kim, J., Shao, W., Schlecht, S.H., Baek, S.Y., Jones, A.K., Ahn, T., Ashton-Miller, J.A., Banaszak Holl, M.M., Wojtys, E.M., 2019. An anterior cruciate ligament failure mechanism. Am. J. Sports Med. 47, 2067–2076. https://doi.org/10.1177/0363546519854450
###### [2] Vila Pouca, M.C.P., Parente, M.P.L., Natal Jorge, R.M., Ashton-Miller, J.A., 2020. Investigating the birth-related caudal maternal pelvic floor muscle injury: the consequences of low cycle fatigue damage. J. Mech. Behav. Biomed. Mater 110, 103956. https://doi.org/10.1016/j.jmbbm.2020.103956.
###### [3] Vila Pouca, M. C. P., et al. “Modeling Permanent Deformation during Low-Cycle Fatigue: Application to the Pelvic Floor Muscles during Labor.” Journal of the Mechanics and Physics of Solids, vol. 164, 2022, p. 104908, https://doi.org/10.1016/j.jmps.2022.104908.

# Example 2 - Rectangular Specimen under biaxial cyclic loading

![biaxial](https://user-images.githubusercontent.com/95075305/170662306-0929fe3b-2ff6-4e65-baf5-0bd01e6e066b.png)


- Only 1/8 of the specimen was modeled and simetry boundary conditions were applied
- The specimen was subjected to 20 loading cycles where each cycle consists of:
    - stretch to 1.2
    - unstretch to 1.1
- In the last cycle, the specimen is allowed to unload without an imposed displacement so that we can see the permanent deformation.

![biaxial_compressed](https://user-images.githubusercontent.com/95075305/170665036-a492ef4f-63fd-4421-a068-50c1d4811fa6.gif)
