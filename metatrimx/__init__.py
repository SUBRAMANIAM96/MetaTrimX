# metatrimx/__init__.py

# Print message when package is imported
print("MetaTrimX package initialized!")

# Define the package version
__version__ = "1.0.0"

# Ensure the main script can be accessed from the package
from . import cli  # If you have a cli.py file inside metatrimx
from . import core  # If you have a core.py file inside metatrimx
