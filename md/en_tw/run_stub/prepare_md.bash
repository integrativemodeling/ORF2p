#!/bin/bash
set -e -v

GMX=/home/arthur/progs/gmx-2022.3-single-openmp-gpu/bin/gmx

echo 1 | ${GMX} pdb2gmx -ff amber14sb_OL15_cufix_zn -f ${1} -o conf.gro -ignh -water tip3p
echo System | ${GMX} editconf -f conf.gro -princ -d 1.5 -o conf_c.gro
${GMX} solvate -cp conf_c.gro -cs spc216.gro -p topol.top -o conf_s.gro
${GMX} grompp -f ions.mdp -c conf_s.gro -p topol.top -o ions.tpr -maxwarn 1
echo SOL | ${GMX} genion -s ions.tpr -o conf_si.gro -p topol.top -neutral -conc 0.15
${GMX} grompp -f minim.mdp -c conf_si.gro -p topol.top -o em.tpr -r conf_si.gro
${GMX} mdrun -v -deffnm em
${GMX} grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro
${GMX} mdrun -v -deffnm nvt
#${GMX} energy -f nvt.edr -o temperature.xvg
${GMX} grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -r nvt.gro
${GMX} mdrun -v -deffnm npt
#${GMX} energy -f npt.edr -o pressure.xvg
${GMX} grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -r npt.gro
${GMX} mdrun -v -deffnm md 
