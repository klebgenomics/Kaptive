import sys
from pathlib import Path

modules_to_exclude = ["kaptive.cli", "kaptive.client", "kaptive.plotting", "kaptive.bgc", "kaptive._version"]
src_dir = Path("src/kaptive")
docs_api_dir = Path("docs/reference")
docs_api_dir.mkdir(parents=True, exist_ok=True)

nav = {}

for py_file in src_dir.rglob("*.py"):
    rel_path = py_file.relative_to(src_dir)
    module_parts = ["kaptive"] + list(rel_path.with_suffix("").parts)
    if module_parts[-1] == "__init__":
        module_parts.pop()
    module_name = ".".join(module_parts)
    
    if module_name in modules_to_exclude or any(module_name.startswith(ex + ".") for ex in modules_to_exclude):
        continue
        
    # Write the markdown file
    md_file_path = docs_api_dir / rel_path.with_suffix(".md")
    if md_file_path.name == "__init__.md":
        md_file_path = md_file_path.with_name("index.md")
    
    md_file_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(md_file_path, "w") as f:
        f.write(f"# {module_name}\n\n::: {module_name}\n")
        
    print(f"Created {md_file_path} for {module_name}")
