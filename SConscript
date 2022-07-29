#!/bin/env python3
import os

Import("exportEnv")
envClone = exportEnv.Clone()

# source_directories = ['core', 'higgsino', ...]
source_directories = next(os.walk('src'))[1]
if len(source_directories) == 0:
  source_directories = ['./']

envClone.Append(CPPPATH = ['#/inc'])

# Determine names for trees
# tree_variable_filenames = ['pico']
tree_variable_filenames = [os.path.basename(str(x)) for x in Glob('#/txt/variables/*')]
tree_variable_files = []
tree_generated_files = []
tree_generated_files_mark = []
for tree_name in tree_variable_filenames: # tree_name = 'pico'
  tree_variable_files.append('txt/variables/'+tree_name)
  tree_generated_files.extend(['src/'+tree_name+'_tree.cpp',
                               'inc/'+tree_name+'_tree.hpp',
                              ])
tree_generated_files_mark = ['#/'+term for term in tree_generated_files]

# Make tree generator
tree_generator_file = 'src/generate_tree_classes.cxx'
tree_generator_exe = 'kernel/'+envClone['kernel']+'/run/generate_tree_classes.exe'
tree_generator = envClone.Program('#/'+tree_generator_exe, Glob(tree_generator_file))

# Run tree generator
run_tree_generator = envClone.Command(tree_generated_files_mark, [], './'+tree_generator_exe+' '+' '.join(tree_variable_filenames))
# input: tree_generator, tree_variable_file
envClone.Depends(run_tree_generator, tree_generator_file)
envClone.Depends(run_tree_generator, tree_variable_files)
envClone.Depends(run_tree_generator, tree_generator)

# Make binaries for every directory
exclude_files = [tree_generator_file] + tree_generated_files
for source_directory in source_directories: # source_directories = 'core'
  library_objects = []
  # Build baby objects
  for source_file in tree_generated_files:
    if 'hpp' in source_file: continue
    source_object = envClone.Object(source_file.replace('cpp','os'), '#/'+source_file)
    library_objects.append(source_object)
  # Build non-baby core objects
  for lib_file in Glob("src/"+source_directory+"*.cpp", exclude=tree_generated_files): 
    source_object = envClone.Object(lib_file)
    library_objects.append(source_object)
  # Make program
  for source_file in Glob("src/"+source_directory+"/*.cxx", exclude=exclude_files):
    # Build object
    source_object = envClone.Object(source_file)
    # Make binary
    program = envClone.Program('#/kernel/'+envClone['kernel']+'/run/'+source_directory+'/${SOURCE.filebase}.exe', [source_object]+library_objects)
    # Make script that links to binary
    source_basename = os.path.splitext(os.path.basename(str(source_file)))[0]
    source_script_path = 'run/'+source_directory+'/'+source_basename+'.exe'
    source_script_path_mark = "#/"+source_script_path
    envClone.Command(source_script_path_mark, program, './scripts/make_run_scripts.py '+source_script_path)
