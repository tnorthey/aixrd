***,H2

basis=sto-3g;
geometry={
2

H	0.000  0.000  0.000
H	0.740  0.000  0.000
}

hf
put,molden,h2.mld;

