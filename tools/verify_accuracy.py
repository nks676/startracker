import json
import argparse
import numpy as np
from scipy.spatial.transform import Rotation as R

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", default="../data/test_0/truth.json")
    parser.add_argument("--est", nargs=4, type=float, required=True, help="Estimated quat [x, y, z, w]")
    args = parser.parse_args()
    
    with open(args.truth, 'r') as f:
        truth = json.load(f)
        
    q_truth = truth['quaternion_xyzw']
    R_truth = R.from_quat(q_truth)
    R_est = R.from_quat(args.est)
    
    # Normalize our estimate before conversion
    norm = np.linalg.norm(args.est)
    est_q = np.array(args.est) / norm
    R_est = R.from_quat(est_q)
    
    # Calculate angular distance
    # The minimum angle separating the two rotations
    # Since quaternions q and -q represent same rotation, distance handles it:
    R_err = R_truth * R_est.inv()
    angle_err = R_err.magnitude()
    
    print(f"Ground Truth Quat: {q_truth}")
    print(f"Estimated Quat:    {list(est_q)}")
    print(f"Angle error: {np.degrees(angle_err):.4f} degrees")
    if np.degrees(angle_err) < 2.0:
        print("Verification PASSED!")
    else:
        print("Verification FAILED: Error too large.")
        
if __name__ == "__main__":
    main()
