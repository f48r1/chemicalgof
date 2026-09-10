
def cd_parent():
    import os
    from pathlib import Path

    if 'notebook_file_dir' in globals():
        print(f"[INFO] word directory already set to {os.getcwd()}")
        return

    notebook_path_dir = Path(__file__).resolve()

    globals()['notebook_file_dir'] = str(notebook_path_dir.parent)

    new_workdir = notebook_path_dir.parent.parent
    os.chdir(new_workdir)

    print(f"[INFO] Working directory set to: {new_workdir}")
