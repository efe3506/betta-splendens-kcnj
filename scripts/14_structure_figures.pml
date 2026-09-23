# Structure figures (PyMOL):  pymol -cq scripts/14_structure_figures.pml
# All three figures follow the same rules: selectivity filter (TIGYG) as gold sticks, N and C
# termini labelled, same orientation, residues with pLDDT < 70 removed, except in Fig. S4,
# whose subject is the low-confidence termini.
set ray_opaque_background, 1
set cartoon_transparency, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set specular, 0.2
set label_size, 18
set label_position, [2.5, 2.5, 8]
set label_color, black
set label_font_id, 7
set label_outline_color, white
set float_labels, 1
viewport 1500, 1900   # ray asagida 3000x3800 ile cagrilir (600 dpi @ 127 mm)
bg_color white

python
from pymol import cmd

def label_termini(obj, tag):
    """Label the first and last remaining CA atom."""
    resis = []
    cmd.iterate(f'{obj} and name CA', 'resis.append(resi)', space={'resis': resis})
    if not resis:
        return
    cmd.label(f'{obj} and name CA and resi {resis[0]}', f'"N, {tag}"')
    cmd.label(f'{obj} and name CA and resi {resis[-1]}', f'"C, {tag}"')

def label_filter(obj, text):
    cmd.label(f'{obj} and pepseq TIGYG and name CA and resn GLY', f'"{text}"')
    # keep the label on the first glycine only
    resis = []
    cmd.iterate(f'{obj} and pepseq TIGYG and name CA and resn GLY', 'resis.append(resi)',
                space={'resis': resis})
    for r in resis[1:]:
        cmd.label(f'{obj} and name CA and resi {r}', '""')

cmd.extend('label_termini', label_termini)
cmd.extend('label_filter', label_filter)
python end

# ---------------------------------------------------------------
# Fig. 5: kcnj16b (red) superimposed on kcnj16a (grey)
# ---------------------------------------------------------------
load results/06_protein_model/kcnj16a_A0A6P7L4L5.pdb, kcnj16a
load results/06_protein_model/kcnj16b_A0A6P7N863.pdb, kcnj16b
remove kcnj16a and b < 70
remove kcnj16b and b < 70
super kcnj16b, kcnj16a

hide everything
show cartoon
color grey70, kcnj16a
color firebrick, kcnj16b
select filt, (kcnj16a or kcnj16b) and pepseq TIGYG
show sticks, filt and not hydro
color gold, filt
set stick_radius, 0.25
label_termini kcnj16b, kcnj16b
label_filter kcnj16b, TIGYG

orient kcnj16a
turn z, -90
turn y, 15
zoom visible, buffer=4
ray 3000, 3800
png figures/fig6_kcnj16b_overlay.png, dpi=600

# ---------------------------------------------------------------
# Fig. S4: kcnj4b coloured by pLDDT, whole chain
# ---------------------------------------------------------------
delete all
load results/06_protein_model/kcnj4b_A0A6P7NPA5.pdb, kcnj4b
hide everything
show cartoon, kcnj4b
spectrum b, orange_yellow_cyan_blue, kcnj4b, minimum=30, maximum=95
select filt, kcnj4b and pepseq TIGYG
show sticks, filt and not hydro
color gold, filt
label_termini kcnj4b, kcnj4b
label_filter kcnj4b, TIGYG

orient kcnj4b
turn z, -90
turn y, 15
zoom visible, buffer=4
ray 3000, 3800
png figures/figS4_kcnj4b_plddt.png, dpi=600

# ---------------------------------------------------------------
# Fig. S5: kcnj10b (blue) superimposed on kcnj10a (grey)
# ---------------------------------------------------------------
delete all
load results/06_protein_model/kcnj10a_A0A6P7LKZ5.pdb, kcnj10a
load results/06_protein_model/kcnj10b_A0A6P7L185.pdb, kcnj10b
remove kcnj10a and b < 70
remove kcnj10b and b < 70
super kcnj10b, kcnj10a

hide everything
show cartoon
color grey70, kcnj10a
color skyblue, kcnj10b
select filt, (kcnj10a or kcnj10b) and pepseq TIGYG
show sticks, filt and not hydro
color gold, filt
label_termini kcnj10b, kcnj10b
label_filter kcnj10b, TIGYG

orient kcnj10a
turn z, -90
turn y, 15
zoom visible, buffer=4
ray 3000, 3800
png figures/figS5_kcnj10b_overlay.png, dpi=600

print "fig5, figS4, figS5 written to figures/"
