
import os
import sys

# Disable colors for markdown generation
os.environ["NO_COLOR"] = "1"
sys.argv[0] = "kaptive"

from kaptive.cli import Cli, Database, Type

def generate_command_docs(cmd_class, output_path, category):
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

    with open(output_path, 'wt') as doc:
        doc.write(f"""---
title: {cmd_title}
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/terminal
categories:
  - {category}
---

{cmd.description}

{usage.strip()}
""")

def main():
    generate_command_docs(Database, 'docs/db/cli.md', 'Databases')
    generate_command_docs(Type, 'docs/serotyping/cli.md', 'Serotyping')

if __name__ == '__main__':
    main()