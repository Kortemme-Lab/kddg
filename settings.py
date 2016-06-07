import os
import json

from klab.fs.fsio import read_file
from klab.general.structures import NestedBunch

sys_settings = None

def load():
    global sys_settings
    if not sys_settings:
        settings_file = os.path.splitext(os.path.abspath(__file__))[0] + '.json'
        d = json.loads(read_file(settings_file))
        sys_settings = NestedBunch(d)
    return sys_settings

