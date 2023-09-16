#!/usr/bin/env python

# System libraries
import argparse as ag
from pathlib import Path

# IM
# IMPP
import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.mmcif
import IMP.pmi.dof
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.bayesianem.restraint
import IMP.pmi.macros

# Additional IMP/IHM libraries
import ihm.cross_linkers


def run_orf2p(
        output_dir='./output',
        n_frames=10000,
    ):

    print("output_dir: {}".format(output_dir))
    print("n_frames: {}".format(n_frames))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput()
    s.add_protocol_output(po)
    po.system.title = "Modeling LINE1 ORF2p"

    glob_data_dir = '../../data/'
    fasta_file = Path(glob_data_dir, "ORF2p.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    pdb_file = Path(glob_data_dir, "af_apo_final_lowpass_4.5.pdb")
    structures["ORF2p"] = (pdb_file, "A", "orange")

    regions_ = {
        # EN
        'NTD': (1, 7, 'F', 0),
        'EN': (8, 237, 'R', 0),
        # Linker 1
        'LINKER1LOOP0': (238, 249, 'F', 0),
        'LINKER1HELIX': (250, 258, 'R', 0),
        'LINKER1LOOP1': (259, 259, 'F', 0),
        # Tower TW
        'TOWERBASE0': (260, 277, 'R', 1), # density
        'TOWERLOOP0': (278, 283, 'F', 1), # density
        'TOWERBASE1': (284, 310, 'R', 1), # density
        'TOWERLOOP1': (311, 312, 'F', 0),
        # Tip of the tower
        'TOWERTIP': (313, 352, 'R', 0),
        'TOWERHELIX0': (353, 359, 'R', 0),
        'TOWERLOOP2': (360, 361, 'F', 0),
        'TOWERHELIX1': (362, 370, 'R', 0),
        'TOWERLOOP3': (371, 374, 'F', 0),
        'TOWERHELIX2': (375, 381, 'R', 0),
        'TOWERLOOP4': (382, 392, 'F', 0),
        # Fingers + RT
        'RT': (393, 849, 'R', 0),
        # Thumb
        'THUMBLOOP0': (850, 856, 'F', 0),
        'THUMBHELIX0': (857, 862, 'R', 0),
        'THUMBLOOP1': (863, 863, 'F', 0),
        'THUMBHELIX1': (864, 868, 'R', 0),
        'THUMBLOOP2': (869, 872, 'F', 0),
        'THUMB': (873, 955, 'R', 0),
        # Linker2/Wrist
        'LINKER2LOOP0': (956, 959, 'F', 0),
        'WRIST': (960, 1030, 'R', 1), # density
        'LINKER2LOOP1': (1031, 1032, 'F', 0),
        'LINKER2HELIX0': (1033, 1061, 'R', 0),
        # CTD
        'CTDLOOP0': (1062, 1067, 'F', 0),
        'CTD': (1068, 1275, 'R', 0),
    }

    regions = {}
    for k, v in regions_.items():
        regions[k] = (v[0], v[1], v[2], v[3])

    mols = dict()

    pdb_file, chain, color = structures['ORF2p']

    mol = st.create_molecule(
        name='ORF2p',
        sequence=seqs['ORF2p'],
        chain_id='A'
    )
    mols['ORF2p'] = mol

    atom = mol.add_structure(
        pdb_fn=str(pdb_file),
        chain_id=chain,
        soft_check=True
    )

    atom = list(atom)

    for i, r in enumerate(regions.values()):
        b, e, t, d = r

        res = [1]

        if d == 1:
            mol.add_representation(
                residues=atom[b - 1:e],
                setup_particles_as_densities=1,
                resolutions=res,
                color=float(i) / len(regions)
        )
        else:
            mol.add_representation(
                residues=atom[b - 1:e],
            resolutions=res,
            # setup_particles_as_densities=1,
            color=float(i) / len(regions)
            )

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES

    flexible_particles = []

    for k, r in regions.items():
        b, e, t, d = r

        particles_ = IMP.atom.Selection(
            root_hier, residue_indexes=range(b, e + 1)).get_selected_particles()

        if t == 'R':
            dof.create_rigid_body(
                rigid_parts=particles_,
                max_trans = 1.0,
                max_rot = 0.1,
                nonrigid_max_trans = 0.5
                )

        elif t == 'F':
            flexible_particles += particles_
            dof.create_flexible_beads(flexible_particles, max_trans=0.5)

        elif t == 'C':
            pass
        else:
            raise(Exception('Missing region type'))

    output_objects = []

    # Freeze the RT part

    regs_ = [
        'RT', 'THUMB', 'LINKER2HELIX0'
    ]

    for r_ in regs_:
        r = regions[r_]

        b, e, t, d = r

        particles_ = IMP.atom.Selection(
                root_hier, residue_indexes=range(b, e + 1)).get_selected_particles()

        fixed_beads, fixed_rbs = dof.disable_movers(particles_,
                                                [IMP.core.RigidBodyMover,
                                                 IMP.pmi.TransformMover])

    output_objects = []

    # CONNECTIVITY
    for component in mols.keys():
        rbr = IMP.pmi.restraints.stereochemistry.ResidueBondRestraint(
            objects=mols[component],
            distance=3.81,
            strength=100,
        )
        rbr.add_to_model()
        output_objects.append(rbr)

    # EXCLUDED VOLUME
    ev_objects = list()
    for component in mols.keys():
        ev_objects.append(mols[component])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=1,
    )
    ev_r.set_weight(5)
    ev_r.add_to_model()
    output_objects.append(ev_r)

    # Crosslinks
    xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkwc.set_protein1_key("prot1")
    xldbkwc.set_residue1_key("res1")
    xldbkwc.set_protein2_key("prot2")
    xldbkwc.set_residue2_key("res2")
    xldbkwc.set_psi_key("psi")

    # BS3 crosslinks

    BS3_threshold = 27.0
    BS3_resolution = 1

    crosslinks = [
        ## Data from truncated protein
        ('BS3_RT.csv', 'XL_RT'),
        ## Data from both truncated and full-length proteins
        ('BS3_BOTH.csv', 'XL_BOTH'),
        ## Data from full-length protein
        ('BS3_FL.csv', 'XL_FL'),
        ## Thumb/CTD from MD
        ('BS3_MD.csv', 'XL_MD'),
    ]

    for (fname, label) in crosslinks:
        # full path
        fpath = str(Path(glob_data_dir, 'xls', fname))
        xldb_ = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
        xldb_.create_set_from_file(fpath)
        xl_ = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb_,
            length=BS3_threshold,
            resolution=BS3_resolution,
            label=label,
            linker=ihm.cross_linkers.bs3,
            weight=5.0
        )

        xl_.set_psi_is_sampled(True)
        xl_.add_to_model()
        output_objects.append(xl_)
        dof.get_nuisances_from_restraint(xl_)

    # EM
    em = True

    if em:
        sel = IMP.atom.Selection(
            hierarchy=root_hier,
            representation_type=IMP.atom.DENSITIES
        )
        densities = sel.get_selected_particles()
        target_gmm_file = Path(glob_data_dir, "gmm/map_ng55_cut.txt")
        gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(
            densities,
            target_fn=str(target_gmm_file),
            scale_target_to_mass=True,
            weight=5
        )
        gem.add_to_model()
        gem.set_label("B3DEM")
        output_objects.append(gem)

    rex = IMP.pmi.macros.ReplicaExchange(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        number_of_best_scoring_models=10,
        number_of_frames=n_frames,
        test_mode=False
    )

    rex.execute_macro()

    # Finalize the protocol
    po.finalize()
    with open('output.cif', 'w') as fh:
        ihm.dumper.write(fh, [po.system])

if __name__ == "__main__":

    parser = ag.ArgumentParser()
    parser.add_argument('-n', '--nsteps',
        help='Number of steps per replica',
        default=10000, type=int)
    parser.add_argument('-o', '--output_dir',
        help='Path to the output dir',
        default='output', type=Path)
    args = parser.parse_args()

    run_orf2p(
        output_dir=args.output_dir,
        n_frames=args.nsteps,
    )
