from pathlib import Path
import shutil

folder = Path("/mnt/c/Users/gio/Documents/foldseek/version_02/filter/filtered_structures/afdb/mhc_human")
output_folder = Path("/mnt/c/Users/gio/Documents/foldseek/version_02/filter/filtered_structures/afdb/mhc_human_renamed")

output_folder.mkdir(parents=True, exist_ok=True)

for file in folder.glob("*.pdb"):
    if file.is_file():
        name = file.stem
        suffix = file.suffix

        if "_" in name:
            pdb_id, pattern = name.split("-", 1)

            new_pattern = pattern.replace("trimmed_mhc", "mhc_only")
            new_name = f"{pdb_id}_{new_pattern}{suffix}"
            new_path = output_folder / new_name

            # copia com novo nome (sem conflito de case do Windows)
            if not new_path.exists():
                shutil.copy(file, new_path)
                print(f"{file.name} → {new_name}")
            else:
                print(f"Skipped (exists): {new_name}")
