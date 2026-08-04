import os
import sys

from kaptive.cli import Cli

# Disable colors for markdown generation
os.environ["NO_COLOR"] = "1"
sys.argv[0] = "kaptive"


def generate_command_docs(cmd_class, output_path, category, icon) -> None:
    cli = Cli()
    cmd = cmd_class()
    cli.add_command(cmd)

    cmd_title = cmd.name
    if cmd.aliases:
        cmd_title += f" ({', '.join(cmd.aliases)})"

    usage_text = cmd.parser.format_help()
    usage = f"```text\n{usage_text.strip()}\n```\n"

    if cmd.subcommands:
        usage += "\n## Subcommands\n\n"
        for subcmd in cmd.subcommands:
            subcmd_name = subcmd.name
            if subcmd.aliases:
                subcmd_name += f" ({', '.join(subcmd.aliases)})"

            usage += f"### {subcmd_name}\n\n"
            if subcmd.description:
                usage += f"{subcmd.description}\n\n"
            usage += f"```text\n{subcmd.parser.format_help().strip()}\n```\n\n"

    with open(output_path, "w") as doc:
        doc.write(f"""---
title: {cmd_title}
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: {icon}
categories:
  - {category}
---

{cmd.description}

{usage.strip()}
""")


def main() -> None:
    from kaptive.db.cli import Database

    generate_command_docs(Database, "docs/cli/db.md", "Databases", "lucide/database")

    from kaptive.serotyping.cli import Convert, Type

    generate_command_docs(Type, "docs/cli/serotyping.md", "Serotyping", "lucide/syringe")
    generate_command_docs(Convert, "docs/cli/convert.md", "Serotyping", "lucide/arrow-left-right")


if __name__ == "__main__":
    main()
