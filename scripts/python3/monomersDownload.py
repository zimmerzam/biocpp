from ftplib import FTP
import sys
import argparse

parser = argparse.ArgumentParser(description='Download cDNA sequences from "ftp.ensembl.org".')

parser.add_argument('--monomers',nargs='*', default='',
                   help='store the monomers to download (default: all available monomers)')
parser.add_argument('--folder', default='monomers/',
                   help='Folder where to store downloaded files (default: monomers/)')
parser.add_argument('--list', default=False, action='store_true',
                   help='Only print all available monomers (default: False)')
args = parser.parse_args()

if args.folder[-1]!='/':
  args.folder+='/'

ftp = FTP('ftp.wwpdb.org')
ftp.login()
ftp.cwd('pub/pdb/data/monomers/')

monomers = []
if not args.list:
  if len(args.monomers)==0:
    ftp.retrlines('NLST', monomers.append)

  for name in monomers:
    if len(name) > 5:
      continue
    ftp.retrbinary('RETR '+name, open(args.folder+name, 'wb').write)
else:
  monomers = []
  ftp.retrlines('NLST', monomers.append)
  for name in monomers:
    print(name)
ftp.quit()
