Coupling !

---------------------------------------------
New macro.phy reads as:

<phy type_electric="fe2D" type_thermic="fe2D" nature="coupled" res="1e-6" microMSH="micro.msh" microPHY="micro.phy" methodFE2_thermic="2" methodFE2_electric="1">

where:
type_electric and type_thermic = {"dirichlet"; "periodic"; "vonnneumann"; "fe2D"}
nature = {"coupled"; "thermic"; "electric"}

If nature="coupled", all combinations are possible.
Example: fe2D electric + fe2D thermic 
 or dirichlet electric + fe2D thermic,...

If nature="thermic", the code only reads type_thermic.
If nature="electric",only reads type_electric.

!!! CHANGE first line of micro.phy to: <phy type_thermic="periodic">

-----------------------------------------------
