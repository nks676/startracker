import json
import re
import sys
import argparse
import numpy as np
from scipy.spatial.transform import Rotation as R

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", default="../data/test_0/truth.json")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--est", nargs=4, type=float, help="Estimated quat [x, y, z, w]")
    group.add_argument("--result-file", help="Parse quaternion from startracker output file")

    args = parser.parse_args()

    with open(args.truth, 'r') as f:
        truth = json.load(f)

    if args.result_file:
        with open(args.result_file, 'r') as f:
            content = f.read()
        match = re.search(r'Estimated Quaternion: \[(.*?)\]', content)
        if not match:
            print("ERROR: Could not find 'Estimated Quaternion:' in result file")
            sys.exit(1)
        est_q = np.array([float(x.strip()) for x in match.group(1).split(',')])
    else:
        est_q = np.array(args.est)

    q_truth = truth['quaternion_xyzw']
    R_truth = R.from_quat(q_truth)

    est_q = est_q / np.linalg.norm(est_q)
    R_est = R.from_quat(est_q)

    R_err = R_truth * R_est.inv()
    angle_err = R_err.magnitude()

    print(f"Ground Truth Quat: {q_truth}")
    print(f"Estimated Quat:    {list(est_q)}")
    print(f"Angle error: {np.degrees(angle_err):.4f} degrees")
    if np.degrees(angle_err) < 2.0:
        print("Verification PASSED!")
    else:
        print("Verification FAILED: Error too large.")
        sys.exit(1)

if __name__ == "__main__":
    main()
