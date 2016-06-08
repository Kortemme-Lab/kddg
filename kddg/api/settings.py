import os
import json
import sys

from klab import colortext
from klab.fs.fsio import read_file, write_file
from klab.general.structures import NestedBunch

sys_settings = None


def create_template(settings_file):
    write_file(settings_file, json.dumps({
        "database" :
        {
            "username" : "database_username",
            "hostname" : "myserver.mydomain.com",
            "port" : 3306,
            "database" : "database_name",
            "password" : "password_for_database_username",
            "socket" : "path_to_socket_file e.g. /var/lib/mysql/mysql.sock",
            "host_config_name" : "if_my.cnf_is_used"
        },
        "cache" :
        {
            "cache_dir" : "/path/to/cache/files"
        },
        "PDBTM" :
        {
            "xml" : "/path/to/pdbtmall.xml"
        },
        "monomer_api" :
        {
            "prediction_data_path" : "/path/to/monomeric_job_archives"
        },
        "ppi_api" :
        {
            "prediction_data_path" : "/path/to/ppi_job_archives"
        }
    }, sort_keys = True, indent = 4))


def load():
    global sys_settings
    if not sys_settings:
        settings_file = os.path.splitext(os.path.abspath(__file__))[0] + '.json'
        if not os.path.exists(settings_file):
            create_template(settings_file)
            colortext.warning('\nThe settings file {0} needs to be configured. Exiting.\n'.format(settings_file))
            sys.exit(1)
        d = json.loads(read_file(settings_file))
        sys_settings = NestedBunch(d)
    return sys_settings

