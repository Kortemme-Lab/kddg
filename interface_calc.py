#!/bin/python
import os
import sys
import subprocess
import tempfile

rosetta_scripts_xml_file = os.path.join('ddglib', 'score_partners.xml')
output_db3 = 'output.db3'

def rescore_ddg_monomer_pdb(pdb_path, rosetta_scripts_path, chains_to_move, rosetta_database_path = None, scratch_dir = '/tmp'):
    if ',' in chains_to_move:
        chains_to_move = chains_to_move.split(',')
    else:
        chains_to_move = [char for char in chains_to_move]
    tmp_output_dir = tempfile.mkdtemp(prefix='rescore_interface_', dir=scratch_dir)
    command_line = [rosetta_scripts_path]
    if rosetta_database_path:
        command_line.extend( ['-database', rosetta_database_path] )
    command_line.extend( [
        '-parser:protocol',
        os.path.abspath( rosetta_scripts_xml_file ),
        '-s', pdb_path,
        '-inout:dbms:database_name', output_db3,
        '-parser:script_vars', 'chainstomove=%s' % ''.join(chains_to_move),  'currentscorefxn=interface',
    ] )

    out_str = ''
    for arg in command_line:
        out_str += arg + ' '

    rosetta_output = subprocess.check_output(command_line, cwd=tmp_output_dir)#, stderr=subprocess.STDOUT, stdout=subprocess.STDOUT)
    return os.path.join(tmp_output_dir, output_db3)
    
if __name__ == '__main__':
    rescore_ddg_monomer_pdb(sys.argv[1])
