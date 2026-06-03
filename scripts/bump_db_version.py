import argparse
import semver
import tomlkit


def bump_toml_version(file_path: str, bump_type: str):
    # 1. Read the TOML file while preserving all formatting and comments
    with open(file_path, "r") as f:
        doc = tomlkit.load(f)

    if "version" not in doc:
        raise ValueError(f"Could not find 'version' key in {file_path}")

    current_version = str(doc["version"])
    ver = semver.Version.parse(current_version)

    # 2. Calculate the SemVer bump
    if bump_type == 'major':
        new_version = str(ver.bump_major())
    elif bump_type == 'minor':
        new_version = str(ver.bump_minor())
    elif bump_type == 'patch':
        new_version = str(ver.bump_patch())
    else:
        raise ValueError(f"Unknown bump type: {bump_type}")

    # 3. Update the version in the TOML document object
    doc["version"] = new_version

    # 4. Write it safely back to the file
    with open(file_path, "w") as f:
        tomlkit.dump(doc, f)

    # Print ONLY the new version so bash can capture it: VER=$(python bump_version.py ...)
    print(new_version)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to TOML file")
    parser.add_argument("bump", choices=['major', 'minor', 'patch'])
    args = parser.parse_args()

    bump_toml_version(args.file, args.bump)
