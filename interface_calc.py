#!/bin/python
import os
import sys
import subprocess
import tempfile

rosetta_scripts_xml_file = os.path.join('ddglib', 'score_partners.xml')
output_db3 = 'output.db3'

def rescore_ddg_monomer_pdb(pdb_path, rosetta_scripts_path, chains_to_move, rosetta_database_path = None, scratch_dir = '/tmp', score_fxn = 'interface', round_num = None, struct_type = None, extra_flags = []):
    '''Used to rescore ddg monomer pdb on the fly'''
    if ',' in chains_to_move:
        chains_to_move = chains_to_move.split(',')
    else:
        chains_to_move = [char for char in chains_to_move]
    tmp_output_dir = tempfile.mkdtemp(prefix='rescore_interface_', dir=scratch_dir)
    if len(chains_to_move) == 1:
        chains_to_move_str = chains_to_move[0]
    else:
        chains_to_move_str = ','.join(chains_to_move)
    command_line = [rosetta_scripts_path]
    if rosetta_database_path:
        command_line.extend( ['-database', rosetta_database_path] )
    command_line.extend( [
        '-parser:protocol',
        os.path.abspath( rosetta_scripts_xml_file ),
        '-s', pdb_path,
        '-inout:dbms:database_name', output_db3,
        '-parser:script_vars', 'chainstomove=%s' % chains_to_move_str,  'currentscorefxn=%s' % score_fxn,
    ] )
    command_line.extend( extra_flags )

    out_str = ''
    for arg in command_line:
        out_str += arg + ' '

    try:
        rosetta_output = subprocess.check_output(command_line, cwd=tmp_output_dir)#, stderr=subprocess.STDOUT, stdout=subprocess.STDOUT)
    except Exception:
        print 'Rosetta run failed in folder', tmp_output_dir
        print 'Command line:', ' '.join(command_line)
        return (None, None, None)

    if round_num and struct_type:
        return (round_num, struct_type, os.path.join(tmp_output_dir, output_db3))
    else:
        return os.path.join(tmp_output_dir, output_db3)

if __name__ == '__main__':
    rescore_ddg_monomer_pdb(sys.argv[1])
