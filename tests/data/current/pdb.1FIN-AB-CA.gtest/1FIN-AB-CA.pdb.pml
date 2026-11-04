# PyMOL script -- A

# Load the pdb file
load 1FIN-AB-CA.pdb, 1finA
hide all


# ----------------------------------------------------------------------


create MainStruct_1finA, 1finA and chain A 
show cartoon, MainStruct_1finA
color smudge, MainStruct_1finA 
set cartoon_transparency=0.3, MainStruct_1finA



# ----------------------------------------------------------------------


create F146-1finA, 1finA and resi 146 and chain A and not (name c,n,o)
show sticks, F146-1finA
color magenta, F146-1finA
util.cnc F146-1finA


create P171-1finA, 1finA and resi 171 and chain A and not (name c,o)
show sticks, P171-1finA
color magenta, P171-1finA
util.cnc P171-1finA


create L134-1finA, 1finA and resi 134 and chain A and not (name c,n,o)
show sticks, L134-1finA
color magenta, L134-1finA
util.cnc L134-1finA


create W167-1finA, 1finA and resi 167 and chain A and not (name c,n,o)
show sticks, W167-1finA
color magenta, W167-1finA
util.cnc W167-1finA


create I52-1finA, 1finA and resi 52 and chain A and not (name c,n,o)
show sticks, I52-1finA
color magenta, I52-1finA
util.cnc I52-1finA


create E12-1finA, 1finA and resi 12 and chain A and not (name c,n,o)
show sticks, E12-1finA
color magenta, E12-1finA
util.cnc E12-1finA


create L166-1finA, 1finA and resi 166 and chain A and not (name c,n,o)
show sticks, L166-1finA
color magenta, L166-1finA
util.cnc L166-1finA


create P155-1finA, 1finA and resi 155 and chain A and not (name c,o)
show sticks, P155-1finA
color magenta, P155-1finA
util.cnc P155-1finA


create L124-1finA, 1finA and resi 124 and chain A and not (name c,n,o)
show sticks, L124-1finA
color magenta, L124-1finA
util.cnc L124-1finA


create E162-1finA, 1finA and resi 162 and chain A and not (name c,n,o)
show sticks, E162-1finA
color magenta, E162-1finA
util.cnc E162-1finA


create F80-1finA, 1finA and resi 80 and chain A and not (name c,n,o)
show sticks, F80-1finA
color magenta, F80-1finA
util.cnc F80-1finA


create D206-1finA, 1finA and resi 206 and chain A and not (name c,n,o)
show sticks, D206-1finA
color magenta, D206-1finA
util.cnc D206-1finA


create I209-1finA, 1finA and resi 209 and chain A and not (name c,n,o)
show sticks, I209-1finA
color magenta, I209-1finA
util.cnc I209-1finA


create L78-1finA, 1finA and resi 78 and chain A and not (name c,n,o)
show sticks, L78-1finA
color magenta, L78-1finA
util.cnc L78-1finA


create R157-1finA, 1finA and resi 157 and chain A and not (name c,n,o)
show sticks, R157-1finA
color magenta, R157-1finA
util.cnc R157-1finA


create V164-1finA, 1finA and resi 164 and chain A and not (name c,n,o)
show sticks, V164-1finA
color magenta, V164-1finA
util.cnc V164-1finA


create F216-1finA, 1finA and resi 216 and chain A and not (name c,n,o)
show sticks, F216-1finA
color magenta, F216-1finA
util.cnc F216-1finA


create A149-1finA, 1finA and resi 149 and chain A and not (name c,n,o)
show sticks, A149-1finA
color magenta, A149-1finA
util.cnc A149-1finA


create V79-1finA, 1finA and resi 79 and chain A and not (name c,n,o)
show sticks, V79-1finA
color magenta, V79-1finA
util.cnc V79-1finA


create A183-1finA, 1finA and resi 183 and chain A and not (name c,n,o)
show sticks, A183-1finA
color magenta, A183-1finA
util.cnc A183-1finA


create L37-1finA, 1finA and resi 37 and chain A and not (name c,n,o)
show sticks, L37-1finA
color magenta, L37-1finA
util.cnc L37-1finA


create K237-1finA, 1finA and resi 237 and chain A and not (name c,n,o)
show sticks, K237-1finA
color magenta, K237-1finA
util.cnc K237-1finA


create V30-1finA, 1finA and resi 30 and chain A and not (name c,n,o)
show sticks, V30-1finA
color magenta, V30-1finA
util.cnc V30-1finA


create I35-1finA, 1finA and resi 35 and chain A and not (name c,n,o)
show sticks, I35-1finA
color magenta, I35-1finA
util.cnc I35-1finA


create K34-1finA, 1finA and resi 34 and chain A and not (name c,n,o)
show sticks, K34-1finA
color magenta, K34-1finA
util.cnc K34-1finA


create F285-1finA, 1finA and resi 285 and chain A and not (name c,n,o)
show sticks, F285-1finA
color magenta, F285-1finA
util.cnc F285-1finA


create G153-1finA, 1finA and resi 153 and chain A and not (name c,n,o)
show spheres, G153-1finA
color magenta, G153-1finA
util.cnc G153-1finA


create Q211-1finA, 1finA and resi 211 and chain A and not (name c,n,o)
show sticks, Q211-1finA
color magenta, Q211-1finA
util.cnc Q211-1finA


create P45-1finA, 1finA and resi 45 and chain A and not (name c,o)
show sticks, P45-1finA
color magenta, P45-1finA
util.cnc P45-1finA


create R122-1finA, 1finA and resi 122 and chain A and not (name c,n,o)
show sticks, R122-1finA
color magenta, R122-1finA
util.cnc R122-1finA


create G190-1finA, 1finA and resi 190 and chain A and not (name c,n,o)
show spheres, G190-1finA
color magenta, G190-1finA
util.cnc G190-1finA


create K75-1finA, 1finA and resi 75 and chain A and not (name c,n,o)
show sticks, K75-1finA
color magenta, K75-1finA
util.cnc K75-1finA


create L58-1finA, 1finA and resi 58 and chain A and not (name c,n,o)
show sticks, L58-1finA
color magenta, L58-1finA
util.cnc L58-1finA


create E40-1finA, 1finA and resi 40 and chain A and not (name c,n,o)
show sticks, E40-1finA
color magenta, E40-1finA
util.cnc E40-1finA


create V17-1finA, 1finA and resi 17 and chain A and not (name c,n,o)
show sticks, V17-1finA
color magenta, V17-1finA
util.cnc V17-1finA


create A144-1finA, 1finA and resi 144 and chain A and not (name c,n,o)
show sticks, A144-1finA
color magenta, A144-1finA
util.cnc A144-1finA


create P234-1finA, 1finA and resi 234 and chain A and not (name c,o)
show sticks, P234-1finA
color magenta, P234-1finA
util.cnc P234-1finA


create K178-1finA, 1finA and resi 178 and chain A and not (name c,n,o)
show sticks, K178-1finA
color magenta, K178-1finA
util.cnc K178-1finA


create K6-1finA, 1finA and resi 6 and chain A and not (name c,n,o)
show sticks, K6-1finA
color magenta, K6-1finA
util.cnc K6-1finA


create Y159-1finA, 1finA and resi 159 and chain A and not (name c,n,o)
show sticks, Y159-1finA
color magenta, Y159-1finA
util.cnc Y159-1finA


create P253-1finA, 1finA and resi 253 and chain A and not (name c,o)
show sticks, P253-1finA
color magenta, P253-1finA
util.cnc P253-1finA


create Q5-1finA, 1finA and resi 5 and chain A and not (name c,n,o)
show sticks, Q5-1finA
color magenta, Q5-1finA
util.cnc Q5-1finA


create N59-1finA, 1finA and resi 59 and chain A and not (name c,n,o)
show sticks, N59-1finA
color magenta, N59-1finA
util.cnc N59-1finA


create L189-1finA, 1finA and resi 189 and chain A and not (name c,n,o)
show sticks, L189-1finA
color magenta, L189-1finA
util.cnc L189-1finA


create A282-1finA, 1finA and resi 282 and chain A and not (name c,n,o)
show sticks, A282-1finA
color magenta, A282-1finA
util.cnc A282-1finA


create I104-1finA, 1finA and resi 104 and chain A and not (name c,n,o)
show sticks, I104-1finA
color magenta, I104-1finA
util.cnc I104-1finA


create D235-1finA, 1finA and resi 235 and chain A and not (name c,n,o)
show sticks, D235-1finA
color magenta, D235-1finA
util.cnc D235-1finA


create V123-1finA, 1finA and resi 123 and chain A and not (name c,n,o)
show sticks, V123-1finA
color magenta, V123-1finA
util.cnc V123-1finA


create E57-1finA, 1finA and resi 57 and chain A and not (name c,n,o)
show sticks, E57-1finA
color magenta, E57-1finA
util.cnc E57-1finA


create P228-1finA, 1finA and resi 228 and chain A and not (name c,o)
show sticks, P228-1finA
color magenta, P228-1finA
util.cnc P228-1finA


# Viewing options
bg_color black
reset
