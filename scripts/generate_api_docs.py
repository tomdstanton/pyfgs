import shutil
from pathlib import Path

docs_api_dir = Path("docs/reference")

# Clean existing docs to prevent stale files
if docs_api_dir.exists():
    shutil.rmtree(docs_api_dir)
docs_api_dir.mkdir(parents=True, exist_ok=True)

groups = {
    "Models": ["GeneFinder", "HmmGlobal", "HmmLocal", "Model"],
    "Results": ["Gene", "GeneBatch", "Mutation"],
    "IO": [
        "FastaReader",
        "FastqReader",
        "BedWriter",
        "VcfWriter",
        "Gff3Writer",
        "FnaWriter",
        "FaaWriter",
    ],
}

for group_name, classes in groups.items():
    filename = group_name.lower() + ".md"
    md_file_path = docs_api_dir / filename

    with open(md_file_path, "w") as f:
        f.write("---\nicon: lucide/code\n---\n\n")
        f.write(f"# {group_name}\n\n")

        for cls in classes:
            f.write(f"::: pyfgs.{cls}\n\n")

    print(f"Created {md_file_path} for group {group_name}")
