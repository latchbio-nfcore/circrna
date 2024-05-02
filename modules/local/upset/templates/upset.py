#!/usr/bin/env python3

import pandas as pd
import platform
import upsetplot
import matplotlib
import matplotlib.pyplot as plt
import distutils.version


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

df_tools = pd.DataFrame(
    {
        "tool": "${tools.join(' ')}".split(" "),
        "file": "${beds.join(' ')}".split(" ")
    }
)

tool_files = df_tools.groupby("tool")["file"].apply(lambda x: set(x)).to_dict()
tool_ids = {}

for tool, files in tool_files.items():
    df_tool = pd.concat([pd.read_csv(f, sep="\\t", header=None) for f in files])
    tool_ids[tool] = set(df_tool[3].unique())

dataset = upsetplot.from_contents(tool_ids)

upsetplot.plot(dataset, orientation='horizontal', show_counts=True)
plt.savefig("${meta.id}.png")

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "upsetplot": upsetplot.__version__,
        "matplotlib": matplotlib.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))