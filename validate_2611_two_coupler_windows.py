"""Compatibility entrypoint: both device sets now use terminal-only masks."""
import argparse
from pathlib import Path
from validate_2611_coupler_release import validate, render
if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--render',action='store_true')
    a=p.parse_args();validate(a.root,a.baseline)
    if a.render:render(a.root)
