# CV functions archive

def tica_cv(snapshot, tica_mod, f_scheme, element):
    import openpathsampling as ops
    import mdtraj as md
    
    traj = ops.engines.Trajectory([snapshot]).to_mdtraj()
    traj.remove_solvent(inplace=True)
    
    f_traj = md.compute_contacts(traj, scheme=f_scheme)[0]
    tica_traj = tica_mod.transform(f_traj)[0]
    return tica_traj[element]


def circle(snapshot, tica_1, tica_2, center):
    import math
    return math.sqrt((tica_1(snapshot)-center[0])**2 + (tica_2(snapshot)-center[1])**2)


def dbscan_predict(snapshot, dbscan_mod, tica_1, tica_2):
    '''
    Return the label of a new data point using a dbscan model.
    Iterate over core points and assign a new point to the cluster of the first core point that is within its eps distance. 
    '''
    import numpy as np
    ic1, ic2 = tica_1(snapshot), tica_2(snapshot)
    y_new = -1
    id_a = np.where((np.linalg.norm(dbscan_mod.components_-[ic1, ic2], axis=1) < dbscan_mod.eps) == True)[0]
    if len(id_a)!=0: y_new = dbscan_mod.labels_[dbscan_mod.core_sample_indices_[id_a[0]]]
    return y_new


if __name__=='__main__':
    pass