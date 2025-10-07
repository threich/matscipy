import matscipy
import matscipy.opls
import matscipy.io.opls

parameter_file = 'potoff_opls.in'

s = matscipy.io.opls.read_extended_xyz('struct.extxyz')

parameter_data = matscipy.io.opls.read_parameter_file(parameter_file)
cutoffs, atom_data, bond_data, angle_data, dihed_data = parameter_data

s.set_cutoffs(cutoffs)
s.set_atom_data(atom_data)

assert abs(sum(s.get_charges())) < 0.01  # ensure charge neutrality

s.get_bonds(bond_data)
s.get_angles(angle_data)
s.get_dihedrals(dihed_data)

matscipy.io.opls.write_lammps('struct.lammps', s)
