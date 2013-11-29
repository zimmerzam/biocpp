# This script will find the best radius to be used with 'alpha_exposed.cpp'
# in order to have the best estimation for the list of exposed/non-exposed atoms
# You need to compile 'examples/cgal_support/alpha_exposed.cpp' before using this!!

import glob
import os
import subprocess
import numpy as np

getarea_folder = "../result_GetAreaExposed/"
alphaexposed_executable = "../../examples/cgal_support/alpha_exposed"
getarea_prefix = "getarea_"
alphaexposed_prefix = "alphaexpo_"

# find all pdb files in current folder
pdb_id = [filename[len(getarea_folder):-4] for filename in glob.glob(getarea_folder+"*.pdb")]

# for each pdb file read the corresponding get_area file
getarea = {}
for id in pdb_id:
  getarea_file = open(getarea_folder+getarea_prefix+id+".dat")
  for line in getarea_file:
    token = line.split()
    if not token[0].isdigit():
      continue
    try:
      resid = int(token[3])
      expos = True if float(token[4]) > 0. else False
    except:
      resid = int(token[4])
      expos = True if float(token[5]) > 0. else False
    atname= token[1][:3]
    getarea[(id,resid,atname)] = expos
  getarea_file.close()

# for each pdb file and for each possible water radius compare the output
# from alpha_exposed with that of get_area
for radius in np.arange(1.,18., 0.1):
  error = 0
  ok = 0
  cmd = [alphaexposed_executable, '', '']
  alphexp = {}
  for id in pdb_id:
    cmd[1] = getarea_folder+id+".pdb"
    cmd[2] = str(radius)
    
    logfile = open(alphaexposed_prefix+id+".dat", 'w')
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=logfile)
    p.wait()
    logfile.close()

    alphexp_file = open(alphaexposed_prefix+id+".dat", 'r')
    for line in alphexp_file:
      token = line.split()
      if not token[0].isdigit():
        continue
      resid = int(token[3])
      expos = True if float(token[4]) > 0. else False
      atname= token[1][:3]
      alphexp[(id,resid,atname)] = expos
    alphexp_file.close()
  for key in getarea.iterkeys():
    if key[2][0] not in ['C', 'N', 'S', 'O']:
      continue
    if getarea[key]==True and alphexp.has_key(key):
      ok+=1
    elif getarea[key]==False and not alphexp.has_key(key):
      ok+=1
    else:
      error+=1
  print radius, ok, error, float(error)/float(ok+error)

for id in pdb_id:
  os.remove(alphaexposed_prefix+id+".dat")
