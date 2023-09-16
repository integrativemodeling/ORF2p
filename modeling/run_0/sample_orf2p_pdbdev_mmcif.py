import IMP.pmi.mmcif
from pathlib import Path
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.proteomics
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.bayesianem
import IMP.bayesianem.restraint
import IMP.algebra
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
from IMP.pmi.restraints.basic import ExternalBarrier
import IMP
import IMP.core
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import time
import sys
import math
import pandas as pd
import ihm.cross_linkers


def run_orf2(
        output_dir,
        n_frames,
    ):

    print("output_dir: {}".format(output_dir))
    print("n_frames: {}".format(n_frames))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput()
    s.add_protocol_output(po)
    po.system.title = "Modeling LINE1 ORF2p"
    # Add publication
    # po.system.citations.append(ihm.Citation.from_pubmed_id(25161197))

    glob_data_dir = '../data/'
    fasta_file = Path(glob_data_dir, "ORF2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    pdb_file = Path(glob_data_dir, "apo_final_lowpass_4.5.pdb")
    structures["ORF2"] = (pdb_file, "A", "orange")

    # OFFSET = 237
    OFFSET = 0

    regions_ = {
        # EN
        'NTD': (1, 7, 'F', 'L', 0),
        'EN': (8, 237, 'R', None, 0),
        # Linker 1
        'LINKER1LOOP0': (238, 249, 'F', 'L', 0),
        'LINKER1HELIX': (250, 258, 'R', 'H', 0),
        'LINKER1LOOP1': (259, 259, 'F', 'L', 0),
        # Tower TW
        'TOWERBASE0': (260, 277, 'R', 'H', 1), # density
        'TOWERLOOP0': (278, 283, 'F', 'L', 1), # density
        'TOWERBASE1': (284, 310, 'R', 'H', 1), # density
        'TOWERLOOP1': (311, 312, 'F', 'L', 0),
        # Tip of the tower
        'TOWERTIP': (313, 352, 'R', 'H', 0),
        'TOWERHELIX0': (353, 359, 'R', 'H', 0),
        'TOWERLOOP2': (360, 361, 'F', 'L', 0),
        'TOWERHELIX1': (362, 370, 'R', 'H', 0),
        'TOWERLOOP3': (371, 374, 'F', 'L', 0),
        'TOWERHELIX2': (375, 381, 'R', 'H', 0),
        'TOWERLOOP4': (382, 392, 'F', 'L', 0),
        # Fingers + RT
        'RT': (393, 849, 'R', None, 0),
        # Thumb
        'THUMBLOOP0': (850, 856, 'F', 'L', 0),
        'THUMBHELIX0': (857, 862, 'R', 'H', 0),
        'THUMBLOOP1': (863, 863, 'F', 'L', 0),
        'THUMBHELIX1': (864, 868, 'R', 'H', 0),
        'THUMBLOOP2': (869, 872, 'F', 'L', 0),
        'THUMB': (873, 955, 'R', None, 0),
        # Linker2/Wrist
        'LINKER2LOOP0': (956, 959, 'F', 'L', 0),
        'WRIST': (960, 1030, 'R', 'H', 1), # density
        'LINKER2LOOP1': (1031, 1032, 'F', 'L', 0),
        'LINKER2HELIX0': (1033, 1061, 'R', 'H', 0),
        # CTD
        'CTDLOOP0': (1062, 1067, 'F', 'L', 0),
        'CTD': (1068, 1275, 'R', None, 0),
    }

    regions = {}
    for k, v in regions_.items():
        regions[k] = (v[0] - OFFSET, v[1] - OFFSET, v[2], v[3], v[4])

    mols = dict()

    pdb_file, chain, color = structures['ORF2']

    mol = st.create_molecule(
        name='ORF2',
        sequence=seqs['ORF2'],
        chain_id='A'
    )
    mols['ORF2'] = mol

    atom = mol.add_structure(
        pdb_fn=str(pdb_file),
        chain_id=chain,
        soft_check=True
    )

    atom = list(atom)

    for i, r in enumerate(regions.values()):
        b, e, t, ft, d = r

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
        b, e, t, ft, d = r

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
            if ft == 'H':

                dof.create_rigid_body(
                    rigid_parts=particles_,
                    nonrigid_parts=None,
                )

            elif ft == 'L':
                flexible_particles += particles_
                dof.create_flexible_beads(flexible_particles, max_trans=0.5)

        elif t == 'C':
            pass
        else:
            raise(Exception('Missing region type'))

    output_objects = []

    # Freeze the RT-C part

    regs_ = [
        'RT', 'THUMB', 'LINKER2HELIX0'
    ]

    for r_ in regs_:
        r = regions[r_]

        b, e, t, ft, d = r

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
    #    kappa=0.1
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
    # xldbkwc.set_unique_id_key("uid")
    xldbkwc.set_psi_key("psi")
    # xldbkwc.set_id_score_key("score")

    # BS3 crosslinks

    ## Data from truncated protein

    xldb_ = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xldb_.create_set_from_file(str(Path(glob_data_dir, f'bs3_188d_NaCl_unique.csv')))
    xl_ = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xldb_,
                                   length=27.0,
                                   resolution=1.0,
                                   label=f"XL_RT",
                                   linker=ihm.cross_linkers.bs3,
                                    weight=5.0
    )

    xl_.set_psi_is_sampled(True)
    xl_.add_to_model()
    output_objects.append(xl_)
    dof.get_nuisances_from_restraint(xl_)

    ## Data from both truncated and full-length proteins

    xldb_ = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xldb_.create_set_from_file(str(Path(glob_data_dir, f'bs3_both_unique.csv')))
    xl_ = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xldb_,
                                   length=27.0,
                                   resolution=1.0,
                                   label=f"XL_BOTH",
                                   linker=ihm.cross_linkers.bs3,
                                    weight=5.0
    )

    xl_.set_psi_is_sampled(True)
    xl_.add_to_model()
    output_objects.append(xl_)
    dof.get_nuisances_from_restraint(xl_)

    ## Data from full-length protein

    xldb_ = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xldb_.create_set_from_file(str(Path(glob_data_dir, f'bs3_fl_unique.csv')))
    xl_ = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xldb_,
                                   length=27.0,
                                   resolution=1.0,
                                   label=f"XL_FL",
                                   linker=ihm.cross_linkers.bs3,
                                    weight=5.0
    )

    xl_.set_psi_is_sampled(True)
    xl_.add_to_model()
    output_objects.append(xl_)
    dof.get_nuisances_from_restraint(xl_)

    # C-terminal/thumb from MD

    xldb_ = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xldb_.create_set_from_file(str(Path(glob_data_dir, 'md.csv')))
    xl_ = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xldb_,
                                   length=27.0,
                                   resolution=1.0,
                                   label=f"XL_MD",
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
        test_mode=True
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)
    po.finalize()
    with open('output.cif', 'w') as fh:
        ihm.dumper.write(fh, [po.system])

if __name__ == "__main__":
    run_orf2(
        output_dir=sys.argv[1],
        n_frames=int(sys.argv[2]),
    )
