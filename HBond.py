import numpy as np
import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis.base import AnalysisBase, Results
from MDAnalysis.lib.distances import capped_distance, calc_angles
import os

def main():
    dir_root = "/home/tsingularity/Desktop/MD/SCAN/cc.ct/"
    filename = "traj.lammpstrj"

    for conformer in ['CC']:
        dir_source = os.path.join(dir_root, conformer)
        setting(dir_source, filename)


def setting(dir_source, filename):
    dict_sel = {
        'o_w': 'index 0:125',
        'h': 'index 129:382',
    }

    u = mda.Universe(os.path.join(dir_source, filename), topology_format="LAMMPSDUMP", format="LAMMPSDUMP", dt=0.01)

    hbonds = HBond(
        universe=u,
        donors_sel=dict_sel["o_w"],
        hydrogens_sel=dict_sel["h"],
        acceptors_sel=dict_sel["o_w"],
        d_h_cutoff=1.3,
        d_a_cutoff=3.5,
        h_d_a_angle_cutoff=30,
        update_selections=False)
    hbonds.run(
        start = 0,
        stop = 10,
        verbose = True
    )
    np.savetxt('hbonds.txt', hbonds.hbonds)



class HBond(AnalysisBase):

    def __init__(self, universe,
                 donors_sel=None, hydrogens_sel=None, acceptors_sel=None,
                 d_h_cutoff=1.2, d_a_cutoff=3.0, h_d_a_angle_cutoff=30,
                 update_selections=True):
        self.u = universe
        self._trajectory = self.u.trajectory

        self.donors_sel = donors_sel.strip() if donors_sel is not None else donors_sel
        self.hydrogens_sel = hydrogens_sel.strip() if hydrogens_sel is not None else hydrogens_sel
        self.acceptors_sel = acceptors_sel.strip() if acceptors_sel is not None else acceptors_sel

        self.d_h_cutoff = d_h_cutoff
        self.d_a_cutoff = d_a_cutoff
        self.h_d_a_angle = h_d_a_angle_cutoff
        self.update_selections = update_selections
        self.results = Results()


    def _get_dh_pairs(self):
        hydrogens = self.u.select_atoms(self.hydrogens_sel)
        donors = self.u.select_atoms(self.donors_sel)
        donors_indices, hydrogen_indices = capped_distance(
            reference=donors.positions,
            configuration=hydrogens.positions,
            max_cutoff=self.d_h_cutoff,
            box=self.u.dimensions,
            return_distances=False
        ).T
        donors = donors[donors_indices]
        hydrogens = hydrogens[hydrogen_indices]

        return donors, hydrogens

    def _prepare(self):
        self.results.hbonds = [[], [], [], [], [], []]
        self._acceptors = self.u.select_atoms(self.acceptors_sel,
                                              updating=self.update_selections)
        self._donors, self._hydrogens = self._get_dh_pairs()

    def _single_frame(self):
        box = self._ts.dimensions

        if self.update_selections:
            self._donors, self._hydrogens = self._get_dh_pairs()

        d_a_indices, d_a_distances = capped_distance(
            reference=self._donors.positions,
            configuration=self._acceptors.positions,
            max_cutoff=self.d_a_cutoff,
            min_cutoff=1.0,
            box=box,
            return_distances=True,
        )

        tmp_donors = self._donors[d_a_indices.T[0]]
        tmp_hydrogens = self._hydrogens[d_a_indices.T[0]]
        tmp_acceptors = self._acceptors[d_a_indices.T[1]]

        h_d_a_angles = np.rad2deg(
            calc_angles(
                coords1=tmp_hydrogens.positions,
                coords2=tmp_donors.positions,
                coords3=tmp_acceptors.positions,
                box=box
            )
        )
        hbond_indices = np.where(h_d_a_angles < self.h_d_a_angle)[0]

        hbond_donors = tmp_donors[hbond_indices]
        hbond_hydrogens = tmp_hydrogens[hbond_indices]
        hbond_acceptors = tmp_acceptors[hbond_indices]
        hbond_distances = d_a_distances[hbond_indices]
        hbond_angles = h_d_a_angles[hbond_indices]

        self.results.hbonds[0].extend(np.full_like(hbond_donors,
                                      self._ts.frame))
        self.results.hbonds[1].extend(hbond_donors.indices)
        self.results.hbonds[2].extend(hbond_hydrogens.indices)
        self.results.hbonds[3].extend(hbond_acceptors.indices)
        self.results.hbonds[4].extend(hbond_distances)
        self.results.hbonds[5].extend(hbond_angles)

    def _conclude(self):
        self.results.hbonds = np.asarray(self.results.hbonds).T

    @property
    def hbonds(self):
        return self.results.hbonds


# if __name__ == "__main__":
#     main()
