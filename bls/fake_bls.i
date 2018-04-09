
a
b
c
d
e
(....some variables)



MIND THE SEQUENCE BELOW

1        !(this is rt_scalar_flag(1). This will switch ON the next 3 lines, which will read things about the first scalar. Let us say it is temperature)

temp.dat !(File name of the profile containing file for the first scalar. This name will be pointed to by the variable rt_scalar_profile(1))

10       !(this is the Prandtl number for scalar 1. Will be assigned to rt_pr(1))

0.1      !(this is Richardson number for scalar 1. Will be assigned to rt_ri(1))

1	 !(this is rt_scalar_flag(2). This will switch ON the next 3 lines for properties of second scalar. Let us say it is salt)

salt.dat !(File name of the profile containing file for the second scalar. The name will be pointed to by the variable rt_scalar_flag(2) in the second loop)

100      !(this is the Prandtl number for scalar 2. Will be assigned to rt_pr(2))

0.5      !(this is the Richardson number for scalar 2. Will be assigned to rt_ri(2))


