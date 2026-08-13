import os
import sys

# Disable colors for markdown generation
os.environ["NO_COLOR"] = "1"
sys.argv[0] = "pyfgs"

from pyfgs.cli import get_parser


def generate_cli_docs(
    output_path: str = "docs/cli.md", category: str = "CLI", icon: str = "lucide/terminal"
) -> None:
    parser = get_parser()

    description = parser.description.strip() if parser.description else ""
    usage_text = parser.format_help().strip()

    with open(output_path, "w") as doc:
        doc.write(f"""---
title: CLI Usage
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: {icon}
categories:
  - {category}
---

{description}

```text
{usage_text}
```
""")

    print(f"Created {output_path} for CLI reference")


def main() -> None:
    generate_cli_docs()


if __name__ == "__main__":
    main()
