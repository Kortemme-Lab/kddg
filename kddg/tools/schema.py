import sys

from kddg.api.schema import generate_sqlalchemy_definition, test_schema_against_database_instance
from kddg.api.ppi import get_interface as get_ppi_interface


if __name__ == '__main__':

    generate_sqlalchemy_definition(['PPIDataSetDE', 'PPIDataSetWildTypeDE', 'UserPPAnalysisSetDE'])
    #generate_sqlalchemy_definition([''])
    sys.exit(0)
    ppi_api = get_ppi_interface(read_file(os.path.join('..', 'pw')).strip())
    test_schema_against_database_instance(ppi_api.DDG_db)
