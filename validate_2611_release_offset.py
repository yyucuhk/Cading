"""Compatibility entrypoint; current layout has four sets."""
from validate_2611_four_sets import validate, render
from pathlib import Path
import argparse
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,default=Path(__file__).resolve().parent);p.add_argument("--baseline",type=Path);p.add_argument("--render",action="store_true");a=p.parse_args();validate(a.root,a.baseline)
    if a.render:render(a.root)
