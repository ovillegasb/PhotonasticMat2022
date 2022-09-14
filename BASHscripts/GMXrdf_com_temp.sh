# 0-500ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_com_500.xvg -selrpos res_com -seltype res_com -b 0.0 -e 500.0 -surf no -xvg none

# 500-1000ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_com_1000.xvg -selrpos res_com -seltype res_com -b 500.0 -e 1000.0 -surf no -xvg none

# 1000-1500ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_com_1500.xvg -selrpos res_com -seltype res_com -b 1000.0 -e 1500.0 -surf no -xvg none

# 1500-2000ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_com_2000.xvg -selrpos res_com -seltype res_com -b 1500.0 -e 2000.0 -surf no -xvg none

# 2000-2500ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_com_2500.xvg -selrpos res_com -seltype res_com -b 2000.0 -surf no -xvg none


# 0-500ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_surf_500.xvg -selrpos atom -seltype atom -b 0.0 -e 500.0 -surf mol -xvg none

# 500-1000ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_surf_1000.xvg -selrpos atom -seltype atom -b 500.0 -e 1000.0 -surf mol -xvg none

# 1000-1500ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_surf_1500.xvg -selrpos atom -seltype atom -b 1000.0 -e 1500.0 -surf mol -xvg none

# 1500-2000ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_surf_2000.xvg -selrpos atom -seltype atom -b 1500.0 -e 2000.0 -surf mol -xvg none

# 2000-2500ps
echo -e "2\n 3\n" | gmx rdf -bin 0.1 -f traj_comp.xtc -s confout.gro -o rdf_azo_pol_surf_2500.xvg -selrpos atom -seltype atom -b 2000.0 -surf mol -xvg none

echo -e "2\n" | gmx msd -f traj_comp.xtc -s confout.gro -o msd-azo.xvg