
def cd_parent():
    import os
    from pathlib import Path

    # Se la variabile globale esiste già, non fare nulla
    if 'notebook_file_dir' in globals():
        print(f"[INFO] word directory already set to {os.getcwd()}")
        return

    notebook_path_dir = Path(__file__).resolve()

    # Salva il path della directory corrente del notebook come stringa
    globals()['notebook_file_dir'] = str(notebook_path_dir.parent)

    # Imposta la working directory alla cartella padre
    new_workdir = notebook_path_dir.parent.parent
    os.chdir(new_workdir)

    print(f"[INFO] Working directory set to: {new_workdir}")
