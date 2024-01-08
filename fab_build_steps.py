#!/usr/bin/env python3
from logging import getLogger

from fab.build_config import BuildConfig
from fab.steps.analyse import analyse
from fab.steps.compile_fortran import compile_fortran
from fab.steps.find_source_files import find_source_files
from fab.steps.grab.git import git_checkout
from fab.steps.grab.folder import grab_folder
from fab.steps.link import link_exe
from fab.steps.preprocess import preprocess_fortran

logger = getLogger('fab')

if __name__ == '__main__':

    with BuildConfig(project_label='mckpp-f90') as state:
        #git_checkout(state, src='https://github.com/aosprey/mckpp-f90', revision='main', dst_label='src')
        find_source_files(state)
        preprocess_fortran(state, common_flags=['-DOPENMP'])
        analyse(state, find_programs=True)
        compile_fortran(state)
        link_exe(state)
