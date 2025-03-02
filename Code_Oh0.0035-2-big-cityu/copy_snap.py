import os
import re
import shutil

def find_latest_snapshot(intermediate_dir):
    """
    Find the snapshot file with the largest t value in the intermediate_dir.
    """
    snapshot_files = []
    pattern = re.compile(r"snapshot-(\d+\.\d{4})$")  # Match filenames in the format snapshot-%5.4f

    for filename in os.listdir(intermediate_dir):
        match = pattern.match(filename)
        if match:
            t = float(match.group(1))  # Extract the t value
            snapshot_files.append((t, filename))

    if not snapshot_files:
        return None

    # Sort by t value in descending order and select the file with the largest t
    snapshot_files.sort(reverse=True, key=lambda x: x[0])
    latest_snapshot = snapshot_files[0][1]
    return os.path.join(intermediate_dir, latest_snapshot)

def copy_latest_snapshot_to_dump(subdir):
    """
    Process the intermediate directory in the subfolder and copy the latest snapshot file to dump.
    """
    intermediate_dir = os.path.join(subdir, "intermediate")
    if not os.path.exists(intermediate_dir):
        print(f"No intermediate directory found in {subdir}, skipping.")
        return

    latest_snapshot = find_latest_snapshot(intermediate_dir)
    if not latest_snapshot:
        print(f"No snapshot files found in the intermediate directory of {subdir}, skipping.")
        return

    dump_path = os.path.join(subdir, "dump")
    shutil.copy2(latest_snapshot, dump_path)
    print(f"Copied {latest_snapshot} to {dump_path}")

def main():
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Get all subdirectories in the script's directory
    for entry in os.listdir(script_dir):
        subdir_path = os.path.join(script_dir, entry)
        if os.path.isdir(subdir_path) and entry != "intermediate":  # Exclude the intermediate directory itself
            copy_latest_snapshot_to_dump(subdir_path)

if __name__ == "__main__":
    main()