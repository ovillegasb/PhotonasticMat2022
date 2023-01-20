

# index : a 5 11 12 13
gmx make_ndx -f run.tpr -o dih_ndx.ndx
gmx angle -f traj_comp.xtc -n dih_ndx.ndx -type dihedral -ov angaver.xvg -xvg none

# RDF com-com
gmx rdf -bin 0.05 -f traj_comp.xtc -s run.tpr -o rdf_azo_thf_com.xvg -selrpos res_com -seltype res_com -surf no -xvg none -rmax 1.6 -norm rdf