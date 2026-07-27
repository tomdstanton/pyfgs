import os
import sys

# Disable colors for markdown generation
os.environ["NO_COLOR"] = "1"
sys.argv[0] = "pyfgs"

from pyfgs.cli import get_parser


def main():
    parser = get_parser()

    with open("docs/cli.md", "w") as doc:
        doc.write(f"""---
title: CLI Reference
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/terminal
---

```text
{parser.format_help()}
```
""")


if __name__ == "__main__":
    main()
