"""Compatibility entrypoint: current release-margin and setback audit."""
from validate_2611_release_offset import validate, render
if __name__ == "__main__":
    import argparse
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=Path(__file__).resolve().parent)
    p.add_argument("--baseline",type=Path)
    p.add_argument("--render",action="store_true")
    a=p.parse_args();validate(a.root,a.baseline)
    if a.render:render(a.root)
