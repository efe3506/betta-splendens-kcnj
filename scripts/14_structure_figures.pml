
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set specular, 0.2
set stick_radius, 0.25
bg_color white

load results/06_protein_model/kcnj16a_A0A6P7L4L5.pdb, kcnj16a
load results/06_protein_model/kcnj16b_A0A6P7N863.pdb, kcnj16b
remove (kcnj16a or kcnj16b) and b < 70
super kcnj16b, kcnj16a
hide everything
show cartoon
color grey70, kcnj16a
color firebrick, kcnj16b
select filter_a, kcnj16a and pepseq TIGYG
select filter_b, kcnj16b and pepseq TIGYG
show sticks, (filter_a or filter_b) and not hydro
color black, filter_a
color gold, filter_b
orient kcnj16a
turn z, -90
turn y, 15
zoom visible, buffer=2
ray 1500, 1900
png figures/fig5_kcnj16b_overlay.png, dpi=300

delete all
load results/06_protein_model/kcnj4_A0A6P7NPA5.pdb, kcnj4
hide everything
show cartoon, kcnj4
spectrum b, orange_yellow_cyan_blue, kcnj4, minimum=30, maximum=95
select filter_4, kcnj4 and pepseq TIGYG
show sticks, filter_4 and not hydro
color red, filter_4
orient kcnj4
turn z, -90
zoom visible, buffer=2
ray 1500, 1900
png figures/figS4_kcnj4_plddt.png, dpi=300

delete all
load results/06_protein_model/kcnj10a_A0A6P7LKZ5.pdb, kcnj10a
load results/06_protein_model/kcnj10c_A0A6P7L185.pdb, kcnj10c
remove (kcnj10a or kcnj10c) and b < 70
super kcnj10c, kcnj10a
hide everything
show cartoon
color grey70, kcnj10a
color steelblue, kcnj10c
select filter_10, (kcnj10a or kcnj10c) and pepseq TIGYG
show sticks, filter_10 and not hydro
color gold, filter_10
orient kcnj10a
turn z, -90
turn y, 15
zoom visible, buffer=2
ray 1500, 1900
png figures/figS5_kcnj10c_overlay.png, dpi=300
